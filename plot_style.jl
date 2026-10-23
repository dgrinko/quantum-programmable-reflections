# Shared plot style for the Julia notebooks: quiet chrome, 2px lines, a one-hue blue ramp for
# curves indexed by the dimension d, and direct labels at the curve ends.  `include` this file
# after `using Plots`.  (log_log_holevo_information_plot.ipynb repeats the colours in Python.)

const INK       = "#0b0b0b"    # titles and axis labels
const INK_2     = "#52514e"    # tick labels and direct labels
const MUTED     = "#898781"    # axes and reference lines (e.g. an upper bound)
const GRID      = "#e1e0d9"    # hairline grid
const SURFACE   = "#fcfcfb"    # plot background
const SERIES    = ["#2a78d6", "#eb6834"]   # two-series plots: blue, orange

# Ordinal blue ramp, 10 steps from light to dark.
const BLUE_RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
                   "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

# k colours from the ramp, evenly spaced from its light to its dark end: ramp steps for k ≤ 10,
# interpolated between them for more curves.  Neighbouring steps are hard to tell apart, so
# label the curves directly rather than relying on the colour.
function ramp_colors(k)
    k == 1 && return [BLUE_RAMP[5]]
    k <= length(BLUE_RAMP) && return BLUE_RAMP[round.(Int, range(1, length(BLUE_RAMP); length = k))]
    return [get(cgrad(BLUE_RAMP), t) for t in range(0, 1; length = k)]
end

default(
    size = (640, 420), dpi = 150,
    fontfamily = "sans-serif", titlefontsize = 11, guidefontsize = 10, tickfontsize = 8,
    legendfontsize = 8, titlelocation = :left, tickfontcolor = INK_2,
    foreground_color_text = INK, foreground_color_guide = INK, foreground_color_title = INK,
    foreground_color_axis = MUTED, foreground_color_border = MUTED, tick_direction = :out,
    background_color = SURFACE, foreground_color_legend = nothing, background_color_legend = nothing,
    grid = true, gridcolor = GRID, gridalpha = 1, gridlinewidth = 0.5, gridstyle = :solid,
    framestyle = :axes, linewidth = 2, markersize = 4, markerstrokewidth = 1.5,
    markerstrokecolor = SURFACE, left_margin = 4Plots.mm, bottom_margin = 3Plots.mm,
)

# Ticks for a log axis that spans only a decade or two: 1-2-5 values in [lo, hi], refined with
# 1.5-3-4-7 and then 2.5-3.5-6-8 multiples; if that still gives fewer than three, three
# log-spaced values rounded to two digits.
function log_ticks(lo, hi)
    0 < lo < hi || throw(ArgumentError("log_ticks needs 0 < lo < hi, got ($lo, $hi)"))
    ts = Float64[]
    for mantissas in ((1, 2, 5), (1.5, 3, 4, 7), (2.5, 3.5, 6, 8))
        length(ts) >= 3 && break
        decades = floor(Int, log10(lo)):ceil(Int, log10(hi))
        append!(ts, filter(t -> lo <= t <= hi, [m * 10.0^k for k in decades for m in mantissas]))
    end
    length(ts) < 3 && (ts = round.(exp10.(range(log10(lo), log10(hi); length = 3)); sigdigits = 2))
    ts = unique(sort(ts))
    return ts, [replace(string(round(t; sigdigits = 2)), r"\.0(?=e|$)" => "") for t in ts]
end

# Label a curve at its last point, in secondary ink (text never takes the series colour).
# valign = :top hangs the label below the point, for curves that end just under another one.
end_label!(plt, x, y, text; valign = :vcenter) =
    annotate!(plt, x, y, Plots.text("  " * text, 8, INK_2, :left, valign))

# Save a figure as PDF in the repository's plots/ folder.
save_plot(plt, name) = savefig(plt, joinpath(mkpath(joinpath(@__DIR__, "plots")), name * ".pdf"))
