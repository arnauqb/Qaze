include("constants.jl")

function drawline(x1::Int64, y1::Int64, x2::Int64, y2::Int64)
    x=x1
    y=y1
    dx=abs(x2-x1)
    dy=abs(y2-y1)
    s1=sign(x2-x1)
    s2=sign(y2-y1)
    swap=false
    length = max(dx,dy) + 1
    results = zeros(Int64, length, 2)
    results[1,1] = x1
    results[1,2] = y1
    if(dy>dx)
        temp=dx
        dx=dy
        dy=temp
        swap=true
    end
    p=2*dy-dx
    for i in 2:dx
        while(p>=0)
            p=p-2*dx
            if(swap)
                x+=s1
            else
                y+=s2
            end
        end
        p=p+2*dy;
        if(swap)
            y+=s2
        else
            x+=s1
        end
        results[i, 1] = x
        results[i, 2] = y
    end
    results[length,1] = x2
    results[length,2] = y2
    return results
end
