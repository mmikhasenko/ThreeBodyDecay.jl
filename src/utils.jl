struct minusone end
import Base:^
^(x::minusone, n::Int) = isodd(n) ? -1 : 1
macro x_str(s::String)
    minusone()
end