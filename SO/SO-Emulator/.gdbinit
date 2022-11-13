
set disassembly-flavor intel
set history save on

b so_emul

define dd
disas
end

r

define sr
x/4gx $r9
end


define sidd
si
dd
end

