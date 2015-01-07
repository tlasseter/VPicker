package com.Sts.Actions.Wizards.SeismicAttribute;

import com.Sts.Actions.Wizards.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */
public class StsDefineTypePanel extends JPanel
{
    // boolean hilbertChecked = false;
    JPanel jPanel = new JPanel();
    StsGroupBox hilbertGroupBox = new StsGroupBox("Available Attributes");

    // StsBooleanFieldBean hilbertBean = new StsBooleanFieldBean(this, "hilbertChecked", "Hilbert");
    JCheckBox hilbertCheckbox = new JCheckBox();
    JCheckBox hilbertPhaseCheckbox = new JCheckBox();
    JCheckBox hilbertAmplitudeCheckbox = new JCheckBox();
    JCheckBox hilbertFreqCheckbox = new JCheckBox();
    JCheckBox testVelocityCheckbox = new JCheckBox();
    JCheckBox freqEnhanceCheckbox = new JCheckBox();

    public StsDefineTypePanel(StsWizard wizard, StsWizardStep wizardStep)
    {

        try
        {
            jbInit();
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {

    }

    void jbInit() throws Exception
    {
        this.setLayout(new GridBagLayout());
        this.add(hilbertGroupBox);
        hilbertCheckbox.setText("Hilbert");
        hilbertCheckbox.setSelected(true);
        hilbertPhaseCheckbox.setText("Phase");
        hilbertAmplitudeCheckbox.setText("Amplitude Envelope");
        hilbertFreqCheckbox.setText("Frequency");
        freqEnhanceCheckbox.setText("Freq Enhancement");
        testVelocityCheckbox.setText("Test Velocity");
        hilbertGroupBox.gbc.anchor = hilbertGroupBox.gbc.WEST;
        // hilbertGroupBox.add(hilbertBean);
        hilbertGroupBox.add(hilbertCheckbox);
        hilbertGroupBox.add(hilbertPhaseCheckbox);
        hilbertGroupBox.add(hilbertAmplitudeCheckbox);
        hilbertGroupBox.add(hilbertFreqCheckbox);
        hilbertGroupBox.add(testVelocityCheckbox);
        hilbertGroupBox.add(freqEnhanceCheckbox);
    }

    public ArrayList getAttributes()
    {
        ArrayList list = new ArrayList();
        if (hilbertCheckbox.isSelected())
            list.add(StsSeismicVolumeConstructor.HILBERT);
        if (hilbertPhaseCheckbox.isSelected())
            list.add(StsSeismicVolumeConstructor.HILBERT_PHASE);
        if (hilbertAmplitudeCheckbox.isSelected())
            list.add(StsSeismicVolumeConstructor.HILBERT_AMPLITUDE);
        if (hilbertFreqCheckbox.isSelected())
            list.add(StsSeismicVolumeConstructor.HILBERT_FREQ);
        if (freqEnhanceCheckbox.isSelected())
            list.add(StsSeismicVolumeConstructor.FREQ_ENHANCE);
        return list;
    }

    // public void setHilbertChecked(boolean checked) { hilbertChecked = checked; }
    // public boolean getHilbertChecked() { return hilbertChecked; }
}
