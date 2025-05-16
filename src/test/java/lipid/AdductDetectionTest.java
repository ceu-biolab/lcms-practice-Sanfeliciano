package lipid;

import org.junit.Before;
import org.junit.Test;
import java.util.Set;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

public class AdductDetectionTest {


    @Before
    public void setup() {
    }

    @Test
    public void shouldDetectAdductBasedOnMzDifference() {

        // Given two peaks with ~21.98 Da difference (e.g., [M+H]+ and [M+Na]+)
        Peak mH = new Peak(700.500, 100000.0); // [M+H]+
        Peak mNa = new Peak(722.482, 80000.0);  // [M+Na]+
        Lipid lipid = new Lipid(1, "PC 34:1", "C42H82NO8P", LipidType.PC, 34, 1);

        double annotationMZ = 700.49999d;
        double annotationIntensity = 80000.0;
        double annotationRT = 6.5d;
        Annotation annotation = new Annotation(lipid, annotationMZ, annotationIntensity, annotationRT, Set.of(mH, mNa), Ionization.POSITIVE);


        assertNotNull("[M+H]+ should be detected", annotation.getAdduct());
        assertEquals("Adduct inferred from lowest mz in group", "[M+H]+", annotation.getAdduct());

    }


    @Test
    public void shouldDetectLossOfWaterAdduct() {
        Peak mh = new Peak(700.500, 90000.0);        // [M+H]+
        Peak mhH2O = new Peak(682.4894, 70000.0);     // [M+H–H₂O]+, ~18.0106 Da less

        Lipid lipid = new Lipid(1, "PE 36:2", "C41H78NO8P", LipidType.PE, 36, 2);
        Annotation annotation = new Annotation(lipid, mh.getMz(), mh.getIntensity(), 7.5d, Set.of(mh, mhH2O), Ionization.POSITIVE);


        assertNotNull("[M+H]+ should be detected", annotation.getAdduct());

        assertEquals("Adduct inferred from lowest mz in group", "[M+H]+", annotation.getAdduct());
    }

    @Test
    public void shouldDetectDoublyChargedAdduct() {
        // Assume real M = (700.500 - 1.0073) = 699.4927
        // So [M+2H]2+ = (M + 2.0146) / 2 = 350.7536
        Peak singlyCharged = new Peak(700.500, 100000.0);  // [M+H]+
        Peak doublyCharged = new Peak(350.754, 85000.0);   // [M+2H]2+

        Lipid lipid = new Lipid(3, "TG 54:3", "C57H104O6", LipidType.TG, 54, 3);
        Annotation annotation = new Annotation(lipid, singlyCharged.getMz(), singlyCharged.getIntensity(), 10d, Set.of(singlyCharged, doublyCharged), Ionization.POSITIVE);

        assertNotNull("[M+H]+ should be detected", annotation.getAdduct());

        assertEquals("Adduct inferred from lowest mz in group", "[M+H]+", annotation.getAdduct());
    }

    @Test
    public void shouldDetectNegativeAdductBasedOnMzDifference() {

        // Given two peaks with ~36.97 Da difference (e.g., [M–H]– and [M+Cl]–)
        Peak mH = new Peak(700.500, 100000.0); // [M–H]–
        Peak mCl = new Peak(736.4767, 80000.0);  // [M+Cl]–, con ~35.9767 Da diferencia


        Lipid lipid = new Lipid(4, "PI 38:4", "C47H83O13P", LipidType.PI, 38, 4);

        double annotationMZ = 700.49999d;
        double annotationIntensity = 95000.0;
        double annotationRT = 6.8d;

        Annotation annotation = new Annotation(
                lipid,
                annotationMZ,
                annotationIntensity,
                annotationRT,
                Set.of(mH, mCl),
                Ionization.NEGATIVE
        );

        assertNotNull("[M-H]− should be detected", annotation.getAdduct());
        assertEquals("Adduct inferred from lowest mz in group", "[M-H]−", annotation.getAdduct());
    }

    @Test
    public void shouldDetectAdductFromMultiplePeaks() {
        Peak p1 = new Peak(700.500, 100000.0);     // [M+H]+
        Peak p2 = new Peak(722.489, 80000.0);      // [M+Na]+
        Peak p3 = new Peak(350.753, 85000.0);      // [M+2H]2+
        Peak p4 = new Peak(682.489, 70000.0);      // [M+H–H₂O]+

        Lipid lipid = new Lipid(5, "PC 36:4", "C44H80NO8P", LipidType.PC, 36, 4);

        double rt = 6.0;

        Annotation annotation = new Annotation(lipid, p3.getMz(), p3.getIntensity(), rt, Set.of(p1, p2, p3, p4), Ionization.POSITIVE);

        assertNotNull("Adduct should be detected", annotation.getAdduct());
        assertEquals("[M+2H]2+", annotation.getAdduct());
    }

    @Test
    public void shouldDetectNegativeAdductFromMultiplePeaks() {
        double neutralMass = 699.492724;

        Peak p1 = new Peak(neutralMass - 1.0073, 85000.0);   // [M-H]−
        Peak p2 = new Peak(neutralMass - 1.0073 - 18.0106, 60000.0);   // [M-H-H2O]−
        Peak p3 = new Peak(neutralMass + 34.969, 55000.0);   // [M+Cl]−
        Peak p4 = new Peak(neutralMass + 44.998, 58000.0);   // [M+HCOOH-H]−


        Lipid lipid = new Lipid(6, "PG 34:2", "C40H74O10P", LipidType.PG, 34, 2);
        double annotationMz = p2.getMz();
        double intensity = p2.getIntensity();
        double rt = 5.8;

        Annotation annotation = new Annotation(lipid, annotationMz, intensity, rt, Set.of(p1, p2, p3, p4), Ionization.NEGATIVE);

        assertNotNull("Negative adduct should be detected", annotation.getAdduct());
        assertEquals("[M-H-H2O]−", annotation.getAdduct());
    }
    @Test
    public void shouldDetectDimerAdductAmongNegativeOptions() {

        Peak p1 = new Peak(349.2427, 80000.0);   // [M-H]⁻
        Peak p2 = new Peak(331.2321, 60000.0);   // [M-H-H2O]⁻
        Peak p3 = new Peak(385.2190, 55000.0);   // [M+Cl]⁻
        Peak p4 = new Peak(395.2480, 58000.0);   // [M+HCOOH-H]⁻
        Peak p5 = new Peak(699.4927, 70000.0);   // [2M-H]⁻

        Lipid lipid = new Lipid(7, "FA 18:1", "C18H34O2", LipidType.FA, 18, 1);


        Annotation annotation = new Annotation(lipid, p5.getMz(),p5.getIntensity(),10.0,Set.of(p1, p2, p3, p4, p5),Ionization.NEGATIVE);

        assertEquals("[2M-H]−", annotation.getAdduct());
    }


}