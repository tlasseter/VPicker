package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;

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

public class StsRgbaVVolumePanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;

    private final static int RED = 0;
    private final static int GREEN = 1;
    private final static int BLUE = 2;
    private final static int ALPHA = 3;

    private StsSeismicVolume[] seismicVolumes = null;
    private StsVirtualVolume[] virtualVolumes = null;

    private int nSeismicVolumes = 0;
    private int nVirtualVolumes = 0;

    private int selectedSeismicOneIndex = 0;
    private int selectedSeismicTwoIndex = 0;
    private int selectedSeismicThreeIndex = 0;
    private int selectedSeismicFourIndex = 0;

    ButtonGroup volumeType = new ButtonGroup();
    JLabel greenVolumeLbl = new JLabel();
    JLabel redVolumeLbl = new JLabel();
    JComboBox greenVolume = new JComboBox();
    StsJPanel jPanel5 = StsJPanel.addInsets();
    JComboBox redVolume = new JComboBox();
    JComboBox blueVolume = new JComboBox();
    JLabel blueVolumeLbl = new JLabel();
    JComboBox alphaVolume = new JComboBox();
    JLabel alphaVolumeLbl = new JLabel();
    JRadioButton rgbaRadio = new JRadioButton();
    JRadioButton hsiaRadio = new JRadioButton();
  ButtonGroup typeGroup = new ButtonGroup();
  JPanel jPanel1 = new JPanel();
  GridBagLayout gridBagLayout1 = new GridBagLayout();
  GridBagLayout gridBagLayout3 = new GridBagLayout();
  JPanel jPanel2 = new JPanel();
  GridBagLayout gridBagLayout2 = new GridBagLayout();

    public StsRgbaVVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
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
        alphaVolume.addItem("None");
        for(int v = 0; v < nSeismicVolumes; v++)
        {
            redVolume.addItem(seismicVolumes[v].getName());
            greenVolume.addItem(seismicVolumes[v].getName());
            blueVolume.addItem(seismicVolumes[v].getName());
            alphaVolume.addItem(seismicVolumes[v].getName());
        }
        virtualVolumes = ((StsVirtualVolumeWizard)wizard).getVirtualVolumes();
        nVirtualVolumes = virtualVolumes.length;
        for(int v = 0; v < nVirtualVolumes; v++)
        {
            redVolume.addItem(virtualVolumes[v].getName());
            greenVolume.addItem(virtualVolumes[v].getName());
            blueVolume.addItem(virtualVolumes[v].getName());
            alphaVolume.addItem(virtualVolumes[v].getName());
        }
        if((nSeismicVolumes + nVirtualVolumes) > 0)
        {
            setRGBSeismicVolume(0, RED);
            setRGBSeismicVolume(0, GREEN);
            setRGBSeismicVolume(0, BLUE);
            setRGBSeismicVolume(0, ALPHA);
        }
    }

    // Set the RGB seismic combobox to the current index
    private void setRGBSeismicVolume(int volumeIndex, int RGBA)
    {
        if(RGBA == RED)
            redVolume.setSelectedIndex(volumeIndex);
        else if(RGBA == GREEN)
            greenVolume.setSelectedIndex(volumeIndex);
        else if(RGBA == BLUE)
            greenVolume.setSelectedIndex(volumeIndex);
        else if(RGBA == ALPHA)
            alphaVolume.setSelectedIndex(volumeIndex);
        else
            System.out.println("Unknown RGBA Type.");
    }

    // Get the RGB seismic combobox to the current index
    private int getRGBASeismicVolume(int RGBA)
    {
        if(RGBA == RED)
            return selectedSeismicOneIndex;
        else if(RGBA == GREEN)
            return selectedSeismicTwoIndex;
        else if(RGBA == BLUE)
            return selectedSeismicThreeIndex;
        else if(RGBA == ALPHA)
            return selectedSeismicFourIndex;
        else
            return -1;
    }

    void jbInit() throws Exception
    {
        redVolume.setEnabled(true);
        greenVolume.setEnabled(true);
        blueVolume.setEnabled(true);
        alphaVolume.setEnabled(true);
        rgbaRadio.setHorizontalAlignment(SwingConstants.CENTER);
        rgbaRadio.setSelected(true);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout1);
        hsiaRadio.setHorizontalAlignment(SwingConstants.CENTER);
        greenVolumeLbl.setHorizontalAlignment(SwingConstants.RIGHT);
        redVolumeLbl.setHorizontalAlignment(SwingConstants.RIGHT);
        blueVolumeLbl.setHorizontalAlignment(SwingConstants.RIGHT);
        alphaVolumeLbl.setHorizontalAlignment(SwingConstants.RIGHT);
        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setLayout(gridBagLayout2);
        jPanel1.add(hsiaRadio,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 9, 2, 7), 89, 0));
        jPanel1.add(rgbaRadio,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 2, 2, 0), 93, -4));
        jPanel5.add(jPanel1);
        jPanel5.add(jPanel2);

        rgbaRadio.addActionListener(this);
        hsiaRadio.addActionListener(this);

        jPanel2.add(redVolume,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 12, 0, 3), 175, 4));
        jPanel2.add(redVolumeLbl,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 9, 0, 0), 25, 10));
        jPanel2.add(blueVolumeLbl,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 9, 0, 0), 21, 10));
        jPanel2.add(greenVolumeLbl,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, -1, 0, 0), 23, 10));
        jPanel2.add(alphaVolumeLbl,  new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, -1, 9, 0), 23, 10));
        jPanel2.add(alphaVolume,  new GridBagConstraints(1, 3, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 12, 9, 3), 175, 4));
        jPanel2.add(greenVolume,  new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 12, 0, 3), 175, 4));
        jPanel2.add(blueVolume,  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 12, 0, 3), 175, 4));
        this.setLayout(gridBagLayout3);

        greenVolumeLbl.setText("Green PostStack3d");
        greenVolumeLbl.setFont(new java.awt.Font("Dialog", 3, 11));
        redVolumeLbl.setText("Red PostStack3d");
        redVolumeLbl.setFont(new java.awt.Font("Dialog", 3, 11));
        greenVolume.addActionListener(this);
        greenVolume.addActionListener(this);
        greenVolume.setToolTipText("Select second seismic volume");
        jPanel5.setBorder(BorderFactory.createEtchedBorder());
        redVolume.setToolTipText("Select first seismic volume");
        redVolume.addActionListener(this);
        redVolume.addActionListener(this);
        blueVolume.setToolTipText("Select second seismic volume");
        blueVolume.addActionListener(this);
        blueVolume.addActionListener(this);
        blueVolumeLbl.setFont(new java.awt.Font("Dialog", 3, 11));
        blueVolumeLbl.setText("Blue PostStack3d");
        alphaVolume.addActionListener(this);
        alphaVolume.addActionListener(this);
        alphaVolume.setToolTipText("Select second seismic volume");
        alphaVolumeLbl.setText("Alpha PostStack3d");
        alphaVolumeLbl.setFont(new java.awt.Font("Dialog", 3, 11));
        rgbaRadio.setText("RGBA");
        hsiaRadio.setText("HSIA");

    typeGroup.add(rgbaRadio);
    typeGroup.add(hsiaRadio);
    this.add(jPanel5,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 2, 4, 5), 0, 14));
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        if(source == rgbaRadio)
            activateTypeFields(((StsRgbaVVolume)wizardStep).RGBA);
        else if(source == hsiaRadio)
            activateTypeFields(((StsRgbaVVolume)wizardStep).HSIA);
    }

    private void activateTypeFields(int type)
    {
        if(type == ((StsRgbaVVolume)wizardStep).RGBA)
        {
            redVolumeLbl.setText("Red PostStack3d");
            greenVolumeLbl.setText("Green PostStack3d");
            blueVolumeLbl.setText("Blue PostStack3d");
            blueVolumeLbl.setText("Alpha PostStack3d");
        }
        else
        {
            redVolumeLbl.setText("Hue PostStack3d");
            greenVolumeLbl.setText("Saturation PostStack3d");
            blueVolumeLbl.setText("Intensity PostStack3d");
            blueVolumeLbl.setText("Alpha PostStack3d");
        }
    }

    public StsSeismicBoundingBox getSelectedSeismicOneVolume()
    {
        if(selectedSeismicOneIndex > seismicVolumes.length-1)
            return virtualVolumes[selectedSeismicOneIndex - seismicVolumes.length];
        else
            return seismicVolumes[selectedSeismicOneIndex];
    }
    public StsSeismicBoundingBox getSelectedSeismicTwoVolume()
    {
        if(selectedSeismicTwoIndex > seismicVolumes.length)
            return virtualVolumes[selectedSeismicTwoIndex - seismicVolumes.length];
        else
            return seismicVolumes[selectedSeismicTwoIndex];
    }
    public StsSeismicBoundingBox getSelectedSeismicThreeVolume()
    {
        if(selectedSeismicThreeIndex > seismicVolumes.length)
            return virtualVolumes[selectedSeismicThreeIndex - seismicVolumes.length];
        else
            return seismicVolumes[selectedSeismicThreeIndex];
    }
    public StsSeismicBoundingBox getSelectedSeismicFourVolume()
    {
        int idx = selectedSeismicFourIndex;
        idx = idx -1;
        if(idx == -1)
            return null;

        if(idx > seismicVolumes.length)
            return virtualVolumes[idx - seismicVolumes.length];
        else
            return seismicVolumes[idx];
    }

    public StsSeismicBoundingBox[] getSelectedVolumes()
    {
        StsSeismicBoundingBox[] volumes = null;
        if(getSelectedSeismicFourVolume() == null)
        {
            volumes = new StsSeismicBoundingBox[3];
            volumes[0] = getSelectedSeismicOneVolume();
            volumes[1] = getSelectedSeismicTwoVolume();
            volumes[2] = getSelectedSeismicThreeVolume();
        }
        else
        {
            volumes = new StsSeismicBoundingBox[4];
            volumes[0] = getSelectedSeismicOneVolume();
            volumes[1] = getSelectedSeismicTwoVolume();
            volumes[2] = getSelectedSeismicThreeVolume();
            volumes[3] = getSelectedSeismicFourVolume();
        }
        return volumes;
    }

    public int getColorType()
    {
        if(rgbaRadio.isSelected())
            return ((StsRgbaVVolume)wizardStep).RGBA;
        else
            return ((StsRgbaVVolume)wizardStep).HSIA;
    }

}
