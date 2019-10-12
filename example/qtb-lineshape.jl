using ThreeBodyDecay

let
    ev = LinRange(2.0, 10.0, 100)
    cal = iRhoQTB.(ev.^2, 3.0, 1.0, 0.2, 1.0)
    # plot(ev, [real.(cal) imag(cal)], lab="")
end

# iRhoQTB(5.0^2, 3, 1, 0.2, 2)
