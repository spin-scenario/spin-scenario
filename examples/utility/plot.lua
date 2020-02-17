output_terminal{type = "qt", font = "Arial,14"} --available types include qt, png, eps and tex.
-- create y vector.
local y ={}
 for i = 1, 20 do
   y[i] = i*i
 end

 -- case 2.
--plot(line(y))

-- create x vector.
local x ={}
 for i = 1, 20 do
   x[i] = i/2
 end

 -- case 3.
plot(line(x, y))

 -- case 4.
local vy = table2vec(y)
plot(line(vy))

-- case 5.
local vx = table2vec(x)
plot(line(vx, vy))

-- case 1.
write("vector.txt", vy)
plot(line("vector.txt"))

plot("title<line plot test> ylabel<amplitude/Hz> legend<mydata> gnuplot<set key outside>", line(y))

plot("title<line plot test> xlabel<time/s> ylabel<amplitude/Hz> legend<mydata> gnuplot<set key outside>", line(x,y,"k lp *"))


local y2 ={}
 for i = 1, 20 do
   y2[i] = (20-i)^2
 end

 -- multi lines plot.
plot("title<line plot test> xlabel<time/s> ylabel<amplitude/Hz> legend<data1;data2> gnuplot<set key center top>", line(x,y,"k lp *"), line(x,y2,"r lp +"))

-- specify a selected color scheme for multi lines.
plot("title<line plot test> xlabel<time/s> ylabel<amplitude/Hz> legend<data1;data2> gnuplot<set key center top> color<Paired>", line(x,y,"lp *"), line(x,y2,"lp +"))


-- 2D map plot for external file.
plot("title<image matrix plot> xlabel<x> ylabel<y> color<Spectral> gnuplot<set palette negative\n set size ratio -1>", map("raw_abs.txt"))

-- 3d plot
plot("title<image matrix plot> xlabel<x> ylabel<y> color<Spectral> gnuplot<set palette negative\n set size ratio -1\n set cbtics 1000>", map("raw_abs.txt", "style<3d>"))

-- 2D map plot for Lua table.
local usr_table ={}
local nrows, ncols = 32, 64
for i = 1, nrows do
  for j = 1, ncols do
  usr_table[(i-1)*ncols+j] =  math.random(100)
  end
end

local m = table2mat(usr_table, nrows, ncols)

plot("title<image matrix plot> xlabel<x> ylabel<y> color<Spectral>", map(m))
