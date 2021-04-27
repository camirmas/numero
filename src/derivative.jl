# Calculates the first order numerical derivative based on the forward, 
# backward, or centered difference.
function derivative(fn, x, method="central", step=.01)
    if method == "central"
        result = (fn(x + step) - fn(x - step)) / (2*step)
    elseif method == "forward"
        result = (fn(x + step) - fn(x)) / step
    elseif method == "backward"
        result = (fn(x) - fn(x - step)) / step
    else
        result = (fn(x + step) - fn(x - step)) / step
    end
    
    return result
end
