import datetime


def get_stop_function(target_porosity: float,
                      attempts_per_one: int = None,
                      total_attempts: int = None,
                      time_per_one: float = None,
                      total_time: float = None,
                      ) -> callable:
    """
    returns stop function for PorousBox with state: PorousBox.State argument that returns True,
    if stop conditions are met (i.e. if target porosity is achieved or too much time/attempts were spent)
    arguments are values at which generation will stop
    all time values are in seconds
    :param target_porosity: target porosity < 1
    :param attempts_per_one: attempts to generate one object
    :param total_attempts: total attempts to achieve target porosity
    :param time_per_one: time to generate one object
    :param total_time: total time spent to try to achieve target porosity
    :return:
    """
    return lambda state: \
        (attempts_per_one is not None and state.cur_obj_attempts >= attempts_per_one) or \
        (total_attempts is not None and state.total_attempts >= total_attempts) or \
        (time_per_one is not None and state.cur_obj_time >= time_per_one) or \
        (total_time is not None and state.total_time >= total_time) or \
        (state.cur_porocity >= target_porosity)
    # probably very bad
