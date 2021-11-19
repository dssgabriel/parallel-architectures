set term pngcairo size 1080, 720
set output "bench/plot.png"
set title "N-Body simulation"
set grid
set ylabel "Latency in cycles"
set xlabel "Simulation iteration"
plot "bench/out_c.dat" w l t "C version", "bench/out_sd.dat" w l t "SSE scalar", "bench/out_pd.dat" w l t "SSE packed", "bench/out_soa.dat" w l t "C SoA"
