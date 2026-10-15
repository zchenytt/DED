using CairoMakie
Makie.update_theme!(
    fonts = (
        regular = "Times New Roman",
        bold = "Times New Roman Bold",
        italic = "Times New Roman Italic",
        bold_italic = "Times New Roman Bold Italic",
    )
)
richs(s) = rich(s, font="SimSun")

# ========== 全球风电累计装机数据（GWEC Global Wind Report 2026） ==========
years = [2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025]
onshore_total = [23.8,30.8,38.9,47.0,58.4,73.3,93.2,119.6,157.4,195.3,234.5,277.7,311.8,361.4,421.3,473.6,521.9,566.9,621.0,708.9,774.2,841.3,945.6,1052.9,1206.7]
offshore_total = [0.1,0.3,0.5,0.6,0.7,0.5,0.7,1.2,1.8,2.9,3.9,5.2,6.8,8.5,11.9,14.2,18.7,23.0,29.2,36.1,55.5,64.3,75.2,83.2,92.5]
# ========== 绘图 ==========
fig = Figure(size=(700, 900););

ax = Axis(fig[1,1],
    xlabel = richs("年份/年"),
    xlabelsize=22,
    xlabelvisible=false,
    xticklabelsvisible=false,
    ylabel = rich("全球风电累计装机容量/", rich("GW", font = "Times New Roman"), font = "SimSun"),
    ylabelsize=22,
    xticks = 2001:3:2025,
    xgridvisible=false,
    ygridvisible=false
)

xlims!(ax, 2000, 2026)
barplot!(ax, years, onshore_total; color=:forestgreen, label= richs("陆上风电"))
barplot!(ax, years, offshore_total; offset=onshore_total, color=:teal, label= richs("海上风电"))
axislegend(ax; position=:lt)
# ========== 数据（不变） ==========
years = [2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025]
onshore_new = [6.4,7.1,7.9,8.1,11.4,14.6,20.1,26.5,37.9,38.0,39.5,43.7,34.6,50.1,60.4,52.7,49.0,45.6,54.6,88.4,73.4,68.9,105.7,109.3,155.3]
offshore_new = [0.1,0.2,0.3,0.1,0.1,0.1,0.2,0.5,0.6,1.0,1.0,1.3,1.7,1.7,3.4,2.3,4.5,4.3,6.2,6.9,20.7,8.9,10.8,8.0,9.3]
# ========== 绘图（改用offset堆叠，修复报错） ==========
ax = Axis(fig[2,1],
    xlabel = richs("年份/年"),
    xlabelsize=20,
    xticks=2001:3:2025,
    ylabelsize=22,
    ylabel = rich("全球风电年度新增装机容量/", rich("GW", font = "Times New Roman"), font = "SimSun"),
    xgridvisible=false,
    ygridvisible=false
)
xlims!(ax, 2000, 2026)
barplot!(ax, years, onshore_new; color=:forestgreen, label= richs("陆上风电") )
barplot!(ax, years, offshore_new; offset=onshore_new, color=:teal, label= richs("海上风电") )
axislegend(ax; position=:lt)
fig
