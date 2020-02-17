
 -- Do something to create the raw vector data e.g. a list of experiments that each one will return a desired value.
 local usr_table ={}
 for i = 1, 100 do
   usr_table[i] = i
 end

local v = table2vec(usr_table)

-- Do something else e.g. save it to file or plot the data.
write("vector.txt", v)