package zad1.at429630.ewolucja.config;

import zad1.at429630.ewolucja.instr.Instruction;
import zad1.at429630.ewolucja.instr.Program;
import zad1.at429630.ewolucja.utils.Randomizer;

import java.util.Arrays;

public class Config {
    private Integer ileTur; // int
    private Integer poczIleRobow; // int
    private Program poczProg; // Program
    private Integer poczEnergia; // int
    private Integer ileDajeJedzenia; // int
    private Integer ileRosnieJedzenie; // int
    private Integer kosztTury; // int
    private Double prPowielania; // double
    private Double ulamekEnergiiRodzica; // double
    private Integer limitPowielania; // int
    private Double prUsunieciaInstr; // double
    private Double prDodanieInstr; // double
    private Instruction[] spisInstr; // Instr[]
    private Double prZmianaInstr; // double
    private Integer coIleWypisz; // int
    private Randomizer randomizer;
    private Integer rozmiarPlanszyX;
    private Integer rozmiarPlanszyY;

    void setIleTur(Integer ileTur) {
        this.ileTur = ileTur;
    }

    void setPoczIleRobow(Integer poczIleRobow) {
        this.poczIleRobow = poczIleRobow;
    }

    void setPoczProg(Program poczProg) {
        this.poczProg = poczProg;
    }

    void setPoczEnergia(Integer poczEnergia) {
        this.poczEnergia = poczEnergia;
    }

    void setIleDajeJedzenie(Integer ileDajeJedzenia) {
        this.ileDajeJedzenia = ileDajeJedzenia;
    }

    void setIleRosnieJedzenie(Integer ileRosnieJedzenie) {
        this.ileRosnieJedzenie = ileRosnieJedzenie;
    }

    void setKosztTury(Integer kosztTury) {
        this.kosztTury = kosztTury;
    }

    void setPrPowielania(Double prPowielania) {
        this.prPowielania = prPowielania;
    }

    void setUlamekEnergiiRodzica(Double ulamekEnergiiRodzica) {
        this.ulamekEnergiiRodzica = ulamekEnergiiRodzica;
    }

    void setLimitPowielania(Integer limitPowielania) {
        this.limitPowielania = limitPowielania;
    }

    void setPrUsunieciaInstr(Double prUsunieciaInstr) {
        this.prUsunieciaInstr = prUsunieciaInstr;
    }

    void setPrDodanieInstr(Double prDodanieInstr) {
        this.prDodanieInstr = prDodanieInstr;
    }

    void setSpisInstr(Instruction[] spisInstr) {
        this.spisInstr = spisInstr;
    }

    void setPrZmianaInstr(Double prZmianaInstr) {
        this.prZmianaInstr = prZmianaInstr;
    }

    void setCoIleWypisz(Integer coIleWypisz) {
        this.coIleWypisz = coIleWypisz;
    }

    void setRandomizer(Randomizer randomizer) {
        this.randomizer = randomizer;
    }

    void setRozmiarPlanszyX(int rozmiarPlanszyX) {
        this.rozmiarPlanszyX = rozmiarPlanszyX;
    }

    void setRozmiarPlanszyY(int rozmiarPlanszyY) {
        this.rozmiarPlanszyY = rozmiarPlanszyY;
    }

    public Integer getIleTur() {
        return ileTur;
    }

    public Integer getRozmiarPlanszyX() {
        return rozmiarPlanszyX;
    }

    public Integer getRozmiarPlanszyY() {
        return rozmiarPlanszyY;
    }

    public Integer getPoczIleRobow() {
        return poczIleRobow;
    }

    public Program getPoczProg() {
        return poczProg;
    }

    public Integer getPoczEnergia() {
        return poczEnergia;
    }

    public Integer getIleDajeJedzenia() {
        return ileDajeJedzenia;
    }

    public Integer getIleRosnieJedzenie() {
        return ileRosnieJedzenie;
    }

    public Integer getKosztTury() {
        return kosztTury;
    }

    public Double getPrPowielania() {
        return prPowielania;
    }

    public Double getUlamekEnergiiRodzica() {
        return ulamekEnergiiRodzica;
    }

    public Integer getLimitPowielania() {
        return limitPowielania;
    }

    public Double getPrUsunieciaInstr() {
        return prUsunieciaInstr;
    }

    public Double getPrDodanieInstr() {
        return prDodanieInstr;
    }

    public Instruction[] getSpisInstr() {
        return spisInstr;
    }

    public Double getPrZmianaInstr() {
        return prZmianaInstr;
    }

    public Integer getCoIleWypisz() {
        return coIleWypisz;
    }

    public Randomizer getRandomizer() {
        return randomizer;
    }

    @Override
    public String toString() {
        return "Config{" +
                "ileTur=" + ileTur +
                ", poczIleRobow=" + poczIleRobow +
                ", poczProg=" + poczProg +
                ", poczEnergia=" + poczEnergia +
                ", ileDajeJedzenia=" + ileDajeJedzenia +
                ", ileRosnieJedzenie=" + ileRosnieJedzenie +
                ", kosztTury=" + kosztTury +
                ", prPowielania=" + prPowielania +
                ", ulamekEnergiiRodzica=" + ulamekEnergiiRodzica +
                ", limitPowielania=" + limitPowielania +
                ", prUsunieciaInstr=" + prUsunieciaInstr +
                ", prDodanieInstr=" + prDodanieInstr +
                ", spisInstr=" + Arrays.toString(spisInstr) +
                ", prZmianaInstr=" + prZmianaInstr +
                ", coIleWypisz=" + coIleWypisz +
                '}';
    }
}
