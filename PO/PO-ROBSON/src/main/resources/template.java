
public static void main(String[]args){
    System.out.println(%function_name%);
}

public static Boolean requireBoolean(Boolean b) {
    return b;
}

public static Boolean requireBoolean(Double value) {
    if (value.equals(0.0))
        return false;
    else if (value.equals(1.0))
        return true;

    throw new RuntimeException("Trying to convert " + value + " to Boolean and this is not bool bro.");
}

public static Double requireDouble(Double value) {
    return value;
}

public static Double requireDouble(Boolean b) {
    return b ? 1.0 : 0.0;
}
