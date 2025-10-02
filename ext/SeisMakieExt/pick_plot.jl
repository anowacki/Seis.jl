using Makie

function Seis.pick_axis(ax)

    # Plot crosshairs
    crosshair_x = Makie.Observable([NaN32])
    crosshair_y = Makie.Observable([NaN32])

    vline_pl = Makie.vlines!(ax, crosshair_x; color=:black, linewidth=0.75)
    hline_pl = Makie.hlines!(ax, crosshair_y; color=:black, linewidth=0.75)

    # Vector of clicks returned (immediately) and filled while picking
    clicks = Makie.Observable(typeof((time=0.0, yvalue=0.0))[])
    # Vectors of the same to plot the picks
    times = Makie.Observable(Float32[])
    yvalues = Makie.Observable(Float32[])

    click_pl = Makie.scatter!(ax, times, yvalues; color=:red, marker='+', markersize=15)

    # Update crosshairs
    crosshair_listener = on(Makie.events(ax).mouseposition) do event
        x, y = Makie.mouseposition(ax)
        x0, y0 = ax.finallimits[].origin
        xwidth, ywidth = ax.finallimits[].widths
        if Makie.is_mouseinside(ax) #(x0 <= x <= x0 + xwidth) && (y0 <= y <= y0 + ywidth)
            crosshair_x[][1] = x
            crosshair_y[][1] = y
            Makie.notify(crosshair_x)
            Makie.notify(crosshair_y)
        else
            if !isnan(crosshair_x[][1])
                crosshair_y[][1] = NaN32
                crosshair_x[][1] = NaN32
                Makie.notify(crosshair_x)
                Makie.notify(crosshair_y)
            end
        end
    end

    # Add, remove or clear picks
    pick_listener = on(Makie.events(ax).keyboardbutton; priority=-10) do event
        if !Makie.is_mouseinside(ax)
            return
        end

        # New pick
        if event.key == Makie.Keyboard.a && event.action == Makie.Keyboard.release
            time, yvalue = Makie.mouseposition(ax)
            push!(clicks[], (; time, yvalue))
            Makie.notify(clicks)

            push!(times[], time)
            push!(yvalues[], yvalue)
            Makie.notify(times)
            Makie.notify(yvalues)
            return Makie.Consume(false)
        # Delete last pick
        elseif event.key == Makie.Keyboard.backspace && event.action == Makie.Keyboard.release
            if !isempty(clicks[])
                pop!(clicks[])
                pop!(times[])
                pop!(yvalues[])
                notify(clicks)
                notify(times)
                notify(yvalues)
            end
        # Stop picking
        elseif event.key == Makie.Keyboard.q && event.action == Makie.Keyboard.release
            crosshair_y[][1] = NaN32
            crosshair_x[][1] = NaN32
            Makie.off(crosshair_listener)
            Makie.off(pick_listener)
        end
    end

    return clicks[]
end
