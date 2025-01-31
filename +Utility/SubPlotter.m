function [] = SubPlotter(fig, rows, colmns, total, data, ylabelnames, xlabelname, title)

figure(fig)

for i = 1:rows
subplot(rows,colmns,total)
plot(data)
