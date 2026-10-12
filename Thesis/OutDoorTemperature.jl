const Os = [
    [28,28,27,27,28,29,31,33,35,36,37,38,38,38,37,36,35,34,32,31,30,29,29,28],# Case 1: Typical hot day (peak ~38°C around 15:00)
    [29,28,28,28,29,31,33,35,37,39,40,41,42,42,41,40,38,36,34,33,32,31,30,29],# Case 2: Extremely hot day (peak ~42°C, late afternoon peak)
    [27,27,27,27,28,29,30,32,33,34,35,35,35,35,34,33,32,31,30,29,28,28,28,27],# Case 3: Milder hot day (peak ~35°C, smooth curve)
    [32,32,31,31,32,34,36,38,40,42,43,44,44,44,43,42,41,40,38,37,36,35,34,33],# Case 4: Heatwave night (doesn’t cool much at night, peak 44°C)
    [27,27,27,27,28,29,30,33,35,37,38,39,39,39,38,37,36,34,32,31,30,29,28,27]# Case 5: Cool morning, sharp rise, peak ~39°C
]

using CairoMakie

Makie.update_theme!(
    fonts = (
        regular = "Times New Roman",
        bold = "Times New Roman Bold",
        italic = "Times New Roman Italic",
        bold_italic = "Times New Roman Bold Italic",
    )
)

f = Figure(size=(300,200); backgroundcolor = :transparent, figure_padding = 0.9);
ax = Axis(f[1,1], ylabel= rich("温度/", rich("°C", font="Times New Roman"), font="SimSun") , xlabel=rich("时段/", rich("h", font="Times New Roman"), font="SimSun"), xticks=1:2:24, xgridvisible = false, ygridvisible = false, backgroundcolor = :transparent)
for i = Os
    lines!(ax, 1:24, i)
end
f

save("temperature.svg", f)

