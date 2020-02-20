function script_transform(file)
    local old_file = io.open(file, 'r')
    local script = old_file:read("*all")
    old_file:close()
    -- http://www.lua.org/manual/5.3/manual.html#pdf-string.gsub
    --script = string.gsub(script, "(%w)(#)(%d)", "%1*%3")
    script = string.gsub(script, "(#)(%d)", "*%2")
    script = string.gsub(script, "(])(#)(%d)", "%1*%3") -- in case of seq block comes from lua table.
    script = string.gsub(script, "(})(#)(%d)", "%1*%3")
    script = string.gsub(script, "(%w)(#)", "%1*1") -- if priority number not set, use 1 as default.
    script = string.gsub(script, "(])(#)", "%1*1") -- if priority number not set, use 1 as default. -- in case of seq block comes from lua table.
    --script = string.gsub(script, "(>)(#)", "%1*1")
    script = string.gsub(script, "(%w)(~)", "%1/1")
    --script = string.gsub(script, "(print)", "Print")
    --script = string.gsub(script, "(<)(-?%d*.?%d*)(>)", "%%%2")
    --script = string.gsub(script, "(<%s*)(')(.*)(')(%s*>)", "%%'%3'")
    local new_file = io.open('temp.lua', 'w')
    new_file:write(script)
    new_file:close()
end

function rf_optimizer_tf(par)
    return ssl.cecilia.new(par)
end

function rf_optimizer(par)
    return ssl.rf_optimizer.new(par)
end

function multi_rf_optimizer(par)
    return ssl.multi_rf_optimizer.new(par)
end

-- seq_block module defines
factory = ssl.seq_block_factory.new()
-- ## seq_block_array ## generates multi-seq-block into a table array, which is very easy manipulated via indexing.
function seq_block_array(block_list)
    return factory:create(block_list)
end

-- ## seq_block_scatter ## generates multi-seq-block into individual variables.
function seq_block_scatter(block_list)
    return table.unpack(seq_block_array(block_list))
end

-- ## seq_block ## generates specific seq-block.
function seq_block(par_list)
    local block = table.unpack(seq_block_array { par_list['type'] })
    block:config(par_list)
    return block
end

function observer(par_list)
    local obser = seq_block_scatter { "observer" }
    obser:config(par_list)
    return obser
end

function shapedRF(par_list)
    local rf = seq_block_scatter { "shapedRF" }
    rf:config(par_list)
    return rf
end

function hardRF(par_list)
    local rf = seq_block_scatter { "idealRF" }
    rf:config(par_list)
    return rf
end

function exprGrad(par_list)
    local grad = seq_block_scatter { "analyticalgrad" }
    grad:config(par_list)
    return grad
end

function trapGrad(par_list)
    local grad = seq_block_scatter { "trapezoidgrad" }
    grad:config(par_list)
    return grad
end

function idealGrad(par_list)
    local grad = seq_block_scatter { "idealgrad" }
    grad:config(par_list)
    return grad
end

function delay(par_list)
    local dl = seq_block_scatter { "delay" }
    dl:config(par_list)
    return dl
end

function acq(par_list)
    local acqu = seq_block_scatter { "acquire" }
    acqu:config(par_list)
    return acqu
end

-- this function is only used internally by developer.
function concurrent(par_list)
    local glue = seq_block_scatter { "concurrent" }
    glue:config(par_list)
    return glue
end

-- this function is only used internally by developer.
function serial(par_list)
    local glue = seq_block_scatter { "serial" }
    glue:config(par_list)
    return glue
end

function phantom(file)
    return ssl.phantom.new(file)
end

function spin_system(par_file_or_table)
    return ssl.spin_system.new(par_file_or_table)
end

function line(...)
    return ssl.line.new(...)
end

function lines(...)
    return ssl.line_series.new(...)
end

function map(...)
    return ssl.map.new(...)
end

function xrange(...)
    return ssl.xrange.new(...)
end

function yrange(...)
    return ssl.yrange.new(...)
end

function write(...)
    return ssl.write(...)
end

function plot(...)
    return ssl.plot(...)
end

function Print(...)
    if type(...) == "userdata" then
        return ssl.print(...)
    --else
        --return _G.print(...)
    end
end

function specgram(...)
    return ssl.specgram(...)
end
