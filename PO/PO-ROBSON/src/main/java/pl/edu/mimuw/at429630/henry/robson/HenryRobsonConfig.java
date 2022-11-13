package pl.edu.mimuw.at429630.henry.robson;

import com.squareup.moshi.adapters.PolymorphicJsonAdapterFactory;
import pl.edu.mimuw.at429630.robson.core.InstructionJsonConfig;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;

public class HenryRobsonConfig extends InstructionJsonConfig {
    @Override
    public PolymorphicJsonAdapterFactory<Instruction> getJsonAdapterFactory() {
        return super.getJsonAdapterFactory()
                /* Henry the rob controlling instructions */
                .withSubtype(LeftInstruction.class, "Lewo")
                .withSubtype(RightInstruction.class, "Prawo")
                .withSubtype(ForwardInstruction.class, "Prosto")
                .withSubtype(SmellInstruction.class, "Wachaj")
                .withSubtype(EatInstruction.class, "Jedz");
    }
}
