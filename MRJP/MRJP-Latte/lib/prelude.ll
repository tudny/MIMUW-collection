
; BEGIN: prelude

%struct.String = type<{ i32, i32, i8* }>

; == String constructor/destructor
; Create a new string on heap with the given length and data
declare %struct.String* @_new_str_from_literal(i8* %0, i32 %1)
; Remove the given string at %0 from heap
declare void @_remove_str(%struct.String* %0)

; == String getters
; Get the length of the given string at %0
declare i32 @_get_str_len(%struct.String* %0)
; Get the data of the given string at %0
declare i8* @_get_str_data(%struct.String* %0)

; == String operations
; Create a new String by concatenating the given strings at %0 and %1
declare %struct.String* @_concat_str(%struct.String* %0, %struct.String* %1)
; Check if the given strings at %0 and %1 are equal - return 1 if true, 0 otherwise
declare i32 @_str_eq(%struct.String* %0, %struct.String* %1)
; Check if the given strings at %0 and %1 are not equal - return 1 if true, 0 otherwise
declare i32 @_str_neq(%struct.String* %0, %struct.String* %1)

; == String stdio
; Print the given string at %0
declare void @printString(%struct.String* %0)
; Read a string from stdin and return it
declare %struct.String* @readString()

; == Int stdio
; Print the given integer at %0
declare void @printInt(i32 %0)
; Read an integer from stdin and return it
declare i32 @readInt()

; String mappers
define i1 @_str_eq_i1(%struct.String* %0, %struct.String* %1) {
  %res = call i32 @_str_eq(%struct.String* %0, %struct.String* %1)
  %mapped = icmp eq i32 %res, 1
  ret i1 %mapped
}

define i1 @_str_neq_i1(%struct.String* %0, %struct.String* %1) {
  %res = call i32 @_str_neq(%struct.String* %0, %struct.String* %1)
  %mapped = icmp eq i32 %res, 1
  ret i1 %mapped
}

; Memory utils
declare i8* @_memory_allocate(i32 %0)
declare i8* @_memory_allocate_many(i32 %0, i32 %1)
declare void @_memory_free(i8* %0)
declare i8* @_memory_allocate_many_strings(i32 %0, i8* %1)

declare void @error()

@str.empty = internal constant [0 x i8] c""

; END: prelude


