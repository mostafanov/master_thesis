function [CT, CP] = coaxialrotor(CT_down, CT_up_new_f, cQ_down, cQ_up)
    function [CT_down, CT_up_new_f, cQ_down] = downstreamroto(CT_down, CT_up_new_f, cQ_down)
    end
function [cQ_up] = upstreamrotor(cQ_up)
    end
 
 CT = (CT_down+ CT_up_new_f)
 CP = cQ_down+ cQ_up
end