package com.Sts.Actions.Wizards.MicroseismicCorrelation;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsComputeSetupPanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;

    private StsSeismicVolume[] seismicVolumes = null;
    private StsVirtualVolume[] virtualVolumes = null;

    private int nSeismicVolumes = 0;
    private int nVirtualVolumes = 0;

    private int selectedSeismicOneIndex = 0;

    JComboBox seismicOneCombo = new JComboBox();
    JLabel seismicOneLabel = new JLabel();

    JComboBox operatorComboBox = new JComboBox();
    JLabel jLabel1 = new JLabel();
    JPanel jPanel3 = new JPanel();
    StsCheckbox floatChk = new StsCheckbox("Float Resolution", "Do you want virtual volume computed with float resolution?");
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    transient public static final int AVG = 0;
    transient public static final int MAX = 1;
    transient public static final int MIN = 2;
    transient public static final int VAL = 3;
    transient public static String[] OPERATORS = {"Average", "Maximum", "Minimum", "Value"};

    public StsComputeSetupPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
        //
        // Add all existing seismic volumes
        //
        seismicVolumes = ((StsMicroseismicCorrelationWizard)wizard).getSeismicVolumes();
        nSeismicVolumes = seismicVolumes.length;
        for(int v = 0; v < nSeismicVolumes; v++)
            seismicOneCombo.addItem(seismicVolumes[v].getName());

        virtualVolumes = ((StsMicroseismicCorrelationWizard)wizard).getVirtualVolumes();
        nVirtualVolumes = virtualVolumes.length;
        for(int v = 0; v < nVirtualVolumes; v++)
            seismicOneCombo.addItem(virtualVolumes[v].getName());

        if((nSeismicVolumes + nVirtualVolumes) > 0)
            setSelectedSeismicOneVolume(0);
        //
        // Add operators
        //
        for(int i=0; i < OPERATORS.length; i++)
            operatorComboBox.addItem(OPERATORS[i]);
    }

    // Set the seismic combobox to the current index
    private void setSelectedSeismicOneVolume(int volumeIndex)
    {
        seismicOneCombo.setSelectedIndex(volumeIndex);
    }

    // Get the current seismic index
    public int getSeismicOneVolumeIndex()
    {
        return selectedSeismicOneIndex;
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout2);
        seismicOneLabel.setFont(new Font("Dialog", 3, 11));
        seismicOneLabel.setText("Seismic Volume");
        seismicOneCombo.setToolTipText("Select Seismic Volume");

        operatorComboBox.setMinimumSize(new Dimension(50, 21));
        operatorComboBox.setToolTipText("Select operator to apply");
        jPanel3.setBorder(BorderFactory.createEtchedBorder());
        jPanel3.setLayout(gridBagLayout1);
        jPanel3.add(seismicOneCombo,     new GridBagConstraints(0, 1, 1, 1, 0.5, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 2, 0, 0), 0, 0));
        jPanel3.add(seismicOneLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 0, 0), 0, 0));
        jPanel3.add(operatorComboBox,   new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        this.add(jPanel3, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            , GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 3, 5, 4), 3, -5));

        seismicOneCombo.addActionListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        if(source == seismicOneCombo)
            selectedSeismicOneIndex = seismicOneCombo.getSelectedIndex();
    }

    public StsSeismicVolume getSelectedSeismicOneVolume()
    {
        if(selectedSeismicOneIndex > seismicVolumes.length-1)
            return virtualVolumes[selectedSeismicOneIndex - seismicVolumes.length];
        else
            return seismicVolumes[selectedSeismicOneIndex];
    }

    public StsSeismicVolume getSelectedVolume()
    {
        return getSelectedSeismicOneVolume();
    }

    public int getOperator()
    {
        return operatorComboBox.getSelectedIndex();
    }

}