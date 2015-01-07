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

public class StsBlendVVolumePanel extends JPanel implements ActionListener
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
    
    StsCheckbox floatChk = new StsCheckbox("Float Resolution", "Do you want virtual volume computed with float resolution?");
    
    JPanel jPanel4 = new JPanel();
    JLabel seismicTwoLabel2 = new JLabel();
    JLabel seismicOneLabel2 = new JLabel();
    JComboBox blendTwo = new JComboBox();
    JComboBox blendOne = new JComboBox();
    JTextField blendCondition = new JTextField();
    JLabel jLabel7 = new JLabel();
    JComboBox logicalComboBox = new JComboBox();
  GridBagLayout gridBagLayout1 = new GridBagLayout();
  GridBagLayout gridBagLayout2 = new GridBagLayout();

    public StsBlendVVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
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

    // Set the seismic combobox to the current index
    private void setSelectedBlendOneVolume(int volumeIndex)
    {
        blendOne.setSelectedIndex(volumeIndex);
    }
    // Set the seismic combobox to the current index
    private void setSelectedBlendTwoVolume(int volumeIndex)
    {
        blendTwo.setSelectedIndex(volumeIndex);
    }

    public void initialize()
    {
        model = wizard.getModel();
        //
        // Add all existing seismic volumes
        //
        seismicVolumes = ((StsVirtualVolumeWizard)wizard).getSeismicVolumes();
        nSeismicVolumes = seismicVolumes.length;
        blendTwo.addItem("None");
        for(int v = 0; v < nSeismicVolumes; v++)
        {
            blendOne.addItem(seismicVolumes[v].getName());
            blendTwo.addItem(seismicVolumes[v].getName());
        }
        virtualVolumes = ((StsVirtualVolumeWizard)wizard).getVirtualVolumes();
        nVirtualVolumes = virtualVolumes.length;
        for(int v = 0; v < nVirtualVolumes; v++)
        {
            blendOne.addItem(virtualVolumes[v].getName());
            blendTwo.addItem(virtualVolumes[v].getName());
        }
        if((nSeismicVolumes + nVirtualVolumes) > 0)
        {
            setSelectedBlendOneVolume(0);
            setSelectedBlendTwoVolume(0);
        }
        //
        // Add operators
        //
        for(int i=0; i < StsBlendedVirtualVolume.LOGICALS.length; i++)
            logicalComboBox.addItem(StsBlendedVirtualVolume.LOGICALS[i]);
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

        jPanel4.setLayout(gridBagLayout1);
        jPanel4.setBorder(BorderFactory.createEtchedBorder());
        seismicTwoLabel2.setFont(new java.awt.Font("Dialog", 3, 11));
        seismicTwoLabel2.setText("PostStack3d Two");
        seismicOneLabel2.setFont(new java.awt.Font("Dialog", 3, 11));
        seismicOneLabel2.setText("PostStack3d One");
        blendTwo.setToolTipText("Select second seismic volume");
        blendTwo.addActionListener(this);
        blendOne.setEnabled(true);
    blendOne.setToolTipText("Select first seismic volume");
        blendOne.addActionListener(this);
        blendCondition.setDoubleBuffered(false);
        blendCondition.setToolTipText("Specify condition");
        blendCondition.setText("1.0");
        jLabel7.setFont(new java.awt.Font("Dialog", 1, 11));
        jLabel7.setHorizontalAlignment(SwingConstants.CENTER);
        jLabel7.setText("WHERE");
        logicalComboBox.setToolTipText("Select logical to apply");
    jPanel4.add(seismicOneLabel2,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 9, 0, 0), 83, 0));
        jPanel4.add(blendOne,  new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 9, 0, 0), 124, 4));
    jPanel4.add(jLabel7,  new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 0, 7, 0), 14, -4));
    jPanel4.add(blendTwo,  new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 9), 124, 4));
    jPanel4.add(seismicTwoLabel2,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 9), 84, 0));
    jPanel4.add(logicalComboBox,   new GridBagConstraints(2, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(7, 2, 0, 9), 124, 4));
    jPanel4.add(blendCondition,   new GridBagConstraints(2, 3, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(6, 1, 6, 9), 129, 4));
    jPanel4.add(floatChk,   new GridBagConstraints(1, 3, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(6, 1, 6, 9), 129, 4));    
    this.add(jPanel4,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 7, 6, 5), 7, -1));
        
        floatChk.addActionListener(this);
        blendOne.addActionListener(this);
        blendTwo.addActionListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        if(source == blendOne)
            selectedSeismicOneIndex = blendOne.getSelectedIndex();
        else if(source == blendTwo)
            selectedSeismicTwoIndex = blendTwo.getSelectedIndex();
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
        return logicalComboBox.getSelectedIndex();
    }

    public double getScalar()
    {
        Double scalar = null;

        scalar = new Double(blendCondition.getText());
        return scalar.doubleValue();
    }

}
