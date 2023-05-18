"""
    RunWaitKill(t::Task, time_out::Int)

Schedule a task t1, check it for results every second, schedule another task, t2 for timing, if 
t2 finished before, terminate t1, else, return results t1 and terminate t2. 

"""
function RunWaitKill(t::Task, time_out::Int)
    schedule(t)
    for secondPassed in 1:time_out
        if istaskdone(t) 
            return t.result
        end
        sleep(1)
    end

end


