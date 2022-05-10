function [y] = MyHeaviside(x)

if x>= 1e-15
    y = 1;
else
    y = 0;
end

end

