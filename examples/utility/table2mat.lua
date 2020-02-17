 -- Do something to create the raw matrix data e.g. array experiments that each one will return a desired value.
 local usr_table ={}
 for i = 1, 4 do
   for j = 1, 16 do
   usr_table[(i-1)*16+j] = (i-1)*16+j
   end
 end

local m = table2mat(usr_table, 4, 16)

-- Do something else e.g. save it to file or plot the data.
write("matrix.txt", m)