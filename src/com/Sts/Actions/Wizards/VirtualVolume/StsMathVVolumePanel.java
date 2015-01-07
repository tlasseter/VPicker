package com.Sts.Actions.Wizards.VirtualVolume;

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

public class StsMathVVolumePanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;

    private StsSeismicVolume[] seismicVolumes = null;
    private StsVirtualVolume[] virtualVolumes = null;

    private int nSeismicVolumes = 0;
    private int nVirtualVolumes = 0;

    private int selectedSeismicOneIndex = 0;
    private int selectedSeismicTwoIndex = 0;

    JComboBox seismicOneCombo = new JComboBox();
    JLabel seismicOneLabel = new JLabel();
    JComboBox seismicTwoCombo = new JComboBox();
    JLabel seismicTwoLabel = new JLabel();

    JComboBox operatorComboBox = new JComboBox();
    JLabel jLabel1 = new JLabel();
    JTextField scalarText = new JTextField();
    JPanel jPanel3 = new JPanel();
    StsCheckbox floatChk = new StsCheckbox("Float Resolution", "Do you want virtual volume computed with float resolution?");
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    public StsMathVVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
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
        seismicVolumes = ((StsVirtualVolumeWizard)wizard).getSeismicVolumes();
        nSeismicVolumes = seismicVolumes.length;
        seismicTwoCombo.addItem("None");
        for(int v = 0; v < nSeismicVolumes; v++)
        {
            seismicOneCombo.addItem(seismicVolumes[v].getName());
            seismicTwoCombo.addItem(seismicVolumes[v].getName());
        }
        virtualVolumes = ((StsVirtualVolumeWizard)wizard).getVirtualVolumes();
        nVirtualVolumes = virtualVolumes.length;
        for(int v = 0; v < nVirtualVolumes; v++)
        {
            seismicOneCombo.addItem(virtualVolumes[v].getName());
            seismicTwoCombo.addItem(virtualVolumes[v].getName());
        }
        if((nSeismicVolumes + nVirtualVolumes) > 0)
        {
            setSelectedSeismicOneVolume(0);
            setSelectedSeismicTwoVolume(0);
        }
        //
        // Add operators
        //
        for(int i=0; i < StsMathVirtualVolume.OPERATORS.length; i++)
            operatorComboBox.addItem(StsMathVirtualVolume.OPERATORS[i]);
    }

    // Set the seismic combobox to the current index
    private void setSelectedSeismicOneVolume(int volumeIndex)
    {
        seismicOneCombo.setSelectedIndex(volumeIndex);
    }
    // Set the seismic combobox to the current index
    private void setSelectedSeismicTwoVolume(int volumeIndex)
    {
        seismicTwoCombo.setSelectedIndex(volumeIndex);
    }

    // Get the current seismic index
    public int getSeismicOneVolumeIndex()
    {
        return selectedSeismicOneIndex;
    }
    // Get the current seismic index
    public int getSeismicTwoVolumeIndex()
    {
        return selectedSeismicTwoIndex;
    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout2);
        seismicOneLabel.setFont(new java.awt.Font("Dialog", 3, 11));
        seismicOneLabel.setText("PostStack3d One");
        seismicTwoLabel.setFont(new java.awt.Font("Dialog", 3, 11));
        seismicTwoLabel.setText("PostStack3d Two");

        jLabel1.setFont(new java.awt.Font("Dialog", 1, 11));
        jLabel1.setHorizontalAlignment(SwingConstants.CENTER);
        jLabel1.setText("OR");
        scalarText.setToolTipText("Specify scalar or select second seismic volume");
        scalarText.setText("1.0");
        seismicOneCombo.setToolTipText("Select first seismic volume");
        seismicTwoCombo.setToolTipText("Select optional second seismic volume or specify scalar");
//        seismicOneCombo.setLightWeightPopupEnabled(false);
//        seismicTwoCombo.setLightWeightPopupEnabled(false);
        operatorComboBox.setMinimumSize(new Dimension(50, 21));
        operatorComboBox.setToolTipText("Select operator to apply");
//        operatorComboBox.setLightWeightPopupEnabled(false);
        jPanel3.setBorder(BorderFactory.createEtchedBorder());
        jPanel3.setLayout(gridBagLayout1);
        jPanel3.add(seismicOneCombo,     new GridBagConstraints(0, 1, 1, 1, 0.5, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 2, 0, 0), 124, 4));
        jPanel3.add(seismicOneLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 0, 0), 82, 3));
        jPanel3.add(operatorComboBox,   new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 50, 4));
        jPanel3.add(seismicTwoCombo,     new GridBagConstraints(2, 1, 1, 1, 0.5, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 8), 124, 4));
        jPanel3.add(seismicTwoLabel, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 50), 44, 4));
        jPanel3.add(scalarText, new GridBagConstraints(2, 3, 1, 1, 1.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 5, 8), 129, 4));
        jPanel3.add(jLabel1, new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 43, 0, 61), 44, 5));
        jPanel3.add(floatChk,     new GridBagConstraints(0, 2, 1, 1, 0.5, 0.0
                ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 2, 0, 0), 124, 4));        
        this.add(jPanel3, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            , GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 3, 5, 4), 3, -5));

        floatChk.addActionListener(this);
        seismicOneCombo.addActionListener(this);
        seismicTwoCombo.addActionListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        if(source == seismicOneCombo)
            selectedSeismicOneIndex = seismicOneCombo.getSelectedIndex();
        else if(source == seismicTwoCombo)
            selectedSeismicTwoIndex = seismicTwoCombo.getSelectedIndex();
        else if(source == floatChk)
        	((StsVirtualVolumeWizard)wizard).setFloatVolume(floatChk.isSelected());
    }

    public StsSeismicVolume getSelectedSeismicOneVolume()
    {
        if(selectedSeismicOneIndex > seismicVolumes.length-1)
            return virtualVolumes[selectedSeismicOneIndex - seismicVolumes.length];
        else
            return seismicVolumes[selectedSeismicOneIndex];
    }

    public StsSeismicVolume getSelectedSeismicTwoVolume()
    {
        int idx = selectedSeismicTwoIndex;
        idx = idx -1;
        if(idx == -1)
            return null;
        if(selectedSeismicTwoIndex > seismicVolumes.length)
            return virtualVolumes[idx - seismicVolumes.length];
        else
            return seismicVolumes[idx];
    }

    public StsSeismicVolume[] getSelectedVolumes()
    {
        StsSeismicVolume[] volumes = null;
        if(getSelectedSeismicTwoVolume() == null)
        {
            volumes = new StsSeismicVolume[1];
            volumes[0] = getSelectedSeismicOneVolume();
        }
        else
        {
            volumes = new StsSeismicVolume[2];
            volumes[0] = getSelectedSeismicOneVolume();
            volumes[1] = getSelectedSeismicTwoVolume();
        }
        return volumes;
    }

    public int getOperator()
    {
        return operatorComboBox.getSelectedIndex();
    }

    public double getScalar()
    {
        Double scalar = null;
        scalar = new Double(scalarText.getText());
        return scalar.doubleValue();
    }
}
