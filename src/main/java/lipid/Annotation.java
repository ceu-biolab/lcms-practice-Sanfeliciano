package lipid;

import java.util.*;
import adduct.*;
/**
 * Class to represent the annotation over a lipid
 */
public class Annotation {

    private final Lipid lipid;
    private final double mz;
    private final double intensity;
    private final double rtMin;
    private String adduct;
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;
    private Ionization ionization;
    private static final double PPMTOLERANCE = 10;

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, Ionization ionization) {
        this(lipid, mz, intensity, retentionTime, Collections.emptySet(),ionization);
    }

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param groupedSignals
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, Set<Peak> groupedSignals,Ionization ionization) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.ionization = ionization;
        this.groupedSignals = new TreeSet<>(groupedSignals);
        this.score = 0;
        this.totalScoresApplied = 0;
        detectAdductFromPeaks();

    }

    public Lipid getLipid() {
        return lipid;
    }

    public double getMz() {
        return mz;
    }

    public double getRtMin() {
        return rtMin;
    }

    public String getAdduct() {
        return adduct;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    public double getIntensity() {
        return intensity;
    }

    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }

    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    public Ionization getIonization() {
        return ionization;
    }

    public void setIonization(Ionization ionization) {
        this.ionization = ionization;
    }


    public void addScore(int delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    public double getNormalizedScore() {
        return (double) this.score / this.totalScoresApplied;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation)) return false;
        Annotation that = (Annotation) o;
        return Double.compare(that.mz, mz) == 0 &&
                Double.compare(that.rtMin, rtMin) == 0 &&
                Objects.equals(lipid, that.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format("Annotation(%s, mz=%.4f, RT=%.2f, adduct=%s, intensity=%.1f, score=%d)",
                lipid.getName(), mz, rtMin, adduct, intensity, score);
    }

   public void detectAdductFromPeaks() {

        String finalAdduct = null;

        if (ionization == Ionization.POSITIVE) {
            for (String adduct1 : AdductList.MAPMZPOSITIVEADDUCTS.keySet()) {
                for (String adduct2 : AdductList.MAPMZPOSITIVEADDUCTS.keySet()) {
                    if (adduct1.equals(adduct2)) continue;

                    for (Peak p1 : groupedSignals) {
                        for (Peak p2 : groupedSignals) {
                            if (p1.equals(p2)) continue;

                            double mz1 = p1.getMz();
                            double mz2 = p2.getMz();

                            Double mass1 = Adduct.getMonoisotopicMassFromMZ(mz1, adduct1);
                            Double mass2 = Adduct.getMonoisotopicMassFromMZ(mz2, adduct2);

                            if (mass1 != null && mass2 != null) {
                                int ppmDiffBetweenMasses = Adduct.calculatePPMIncrement(mass1, mass2);

                                if (ppmDiffBetweenMasses <= PPMTOLERANCE) {
                                    int ppmToMz1 = Adduct.calculatePPMIncrement(mz1, this.mz);
                                    int ppmToMz2 = Adduct.calculatePPMIncrement(mz2, this.mz);

                                    if (ppmToMz1 <= PPMTOLERANCE) {
                                        finalAdduct = adduct1;
                                        this.adduct = finalAdduct;
                                        return;
                                    } else if (ppmToMz2 <= PPMTOLERANCE) {
                                        finalAdduct = adduct2;
                                        this.adduct = finalAdduct;
                                        return;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else if (ionization == Ionization.NEGATIVE) {
            for (String adduct1 : AdductList.MAPMZNEGATIVEADDUCTS.keySet()) {
                for (String adduct2 : AdductList.MAPMZNEGATIVEADDUCTS.keySet()) {
                    if (adduct1.equals(adduct2)) continue;

                    for (Peak p1 : groupedSignals) {
                        for (Peak p2 : groupedSignals) {
                            if (p1.equals(p2)) continue;

                            double mz1 = p1.getMz();
                            double mz2 = p2.getMz();

                            Double mass1 = Adduct.getMonoisotopicMassFromMZ(mz1, adduct1);
                            Double mass2 = Adduct.getMonoisotopicMassFromMZ(mz2, adduct2);

                            if (mass1 != null && mass2 != null) {
                                int ppmDiffBetweenMasses = Adduct.calculatePPMIncrement(mass1, mass2);

                                if (ppmDiffBetweenMasses <= PPMTOLERANCE) {
                                    int ppmToMz1 = Adduct.calculatePPMIncrement(mz1, this.mz);
                                    int ppmToMz2 = Adduct.calculatePPMIncrement(mz2, this.mz);

                                    if (ppmToMz1 <= PPMTOLERANCE) {
                                        finalAdduct = adduct1;
                                        this.adduct = finalAdduct;
                                        return;
                                    } else if (ppmToMz2 <= PPMTOLERANCE) {
                                        finalAdduct = adduct2;
                                        this.adduct = finalAdduct;
                                        return;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        this.adduct = finalAdduct;
    }


}



