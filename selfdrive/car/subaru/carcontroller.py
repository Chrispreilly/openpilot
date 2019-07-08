#from common.numpy_fast import clip
from common.realtime import sec_since_boot
from selfdrive.boardd.boardd import can_list_to_can_capnp
from selfdrive.car import apply_std_steer_torque_limits
from selfdrive.car.subaru import subarucan
from selfdrive.car.subaru.values import CAR, DBC
from selfdrive.can.packer import CANPacker


class CarControllerParams():
  def __init__(self, car_fingerprint):

    if car_fingerprint == CAR.IMPREZA:
      self.STEER_MAX = 3500 #2047              # max_steer 4095
      self.STEER_STEP = 2                # how often we update the steer cmd
      self.STEER_DELTA_UP = 50           # torque increase per refresh, 0.8s to max
      self.STEER_DELTA_DOWN = 70         # torque decrease per refresh
      self.STEER_DRIVER_ALLOWANCE = 60   # allowed driver torque before start limiting
      self.STEER_DRIVER_MULTIPLIER = 10   # weight driver torque heavily
      self.STEER_DRIVER_FACTOR = 1     # from dbc
      
    if car_fingerprint == CAR.CROSSTREK:
      self.STEER_MAX = 3500 #2047              # max_steer 4095
      self.STEER_STEP = 2                # how often we update the steer cmd
      self.STEER_DELTA_UP = 50           # torque increase per refresh, 0.8s to max
      self.STEER_DELTA_DOWN = 70         # torque decrease per refresh
      self.STEER_DRIVER_ALLOWANCE = 60   # allowed driver torque before start limiting
      self.STEER_DRIVER_MULTIPLIER = 10   # weight driver torque heavily
      self.STEER_DRIVER_FACTOR = 1     # from dbc

class CarController(object):
  def __init__(self, car_fingerprint):
    self.start_time = sec_since_boot()
    self.lkas_active = False
    self.steer_idx = 0
    self.apply_steer_last = 0
    self.car_fingerprint = car_fingerprint
    self.apply_steer = 0
    self.actuators_steer = 0
    self.es_distance_cnt = -1
    self.es_lkas_cnt = -1

    # Setup detection helper. Routes commands to
    # an appropriate CAN bus number.
    self.params = CarControllerParams(car_fingerprint)
    print(DBC)
    self.packer = CANPacker(DBC[car_fingerprint]['pt'])

  def update(self, enabled, CS, frame, actuators, pcm_cancel_cmd, visual_alert, left_line, right_line):
    """ Controls thread """

    P = self.params

    # Send CAN commands.
    can_sends = []

    ### STEER ###

    if (frame % P.STEER_STEP) == 0:

      final_steer = actuators.steer if enabled else 0.
      apply_steer = int(round(final_steer * P.STEER_MAX))

      # limits due to driver torque

      apply_steer = int(round(apply_steer))
      apply_steer = apply_std_steer_torque_limits(apply_steer, self.apply_steer_last, CS.steer_torque_driver, P)

      lkas_enabled = enabled and not CS.steer_not_allowed

      if not lkas_enabled:
        apply_steer = 0

      can_sends.append(subarucan.create_steering_control(self.packer, CS.CP.carFingerprint, apply_steer, frame, P.STEER_STEP))

      self.apply_steer_last = apply_steer
      # for dashboard
      self.apply_steer = apply_steer
      self.actuators_steer = actuators.steer

    # generate 1Hz op_active msg for global to enable panda es filtering
    if (frame % 100) == 0:
      can_sends.append(subarucan.create_openpilot_active(self.packer))

    if self.es_distance_cnt != CS.es_distance_msg["Counter"]:
      can_sends.append(subarucan.create_es_distance(self.packer, CS.es_distance_msg, pcm_cancel_cmd))
      self.es_distance_cnt = CS.es_distance_msg["Counter"]

    if self.es_lkas_cnt != CS.es_lkas_msg["Counter"]:
      can_sends.append(subarucan.create_es_lkas(self.packer, CS.es_lkas_msg, visual_alert, left_line, right_line))
      self.es_lkas_cnt = CS.es_lkas_msg["Counter"]
 
    return can_sends 
