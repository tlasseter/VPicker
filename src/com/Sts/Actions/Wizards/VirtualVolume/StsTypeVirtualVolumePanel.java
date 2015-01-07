package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
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

public class StsTypeVirtualVolumePanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    StsGroupBox typeBox = new StsGroupBox("Type");
    ButtonGroup typeGrp = new ButtonGroup();

    JRadioButton mathBtn = new JRadioButton("Math");
    JRadioButton blendedBtn = new JRadioButton("Blended");
    JRadioButton xplotBtn = new JRadioButton("Crossplot");
    JRadioButton filterBtn = new JRadioButton("Filter");
    JRadioButton colorBtn = new JRadioButton("Color");
    JRadioButton sensorBtn = new JRadioButton("Sensor");

    StsStringFieldBean nameField = new StsStringFieldBean();
    JTextArea typeDescription = new JTextArea();
    ButtonGroup typeGroup = new ButtonGroup();

    int type = StsVirtualVolume.SEISMIC_MATH;
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    JPanel panel = new JPanel(new GridBagLayout());

    public StsTypeVirtualVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
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
    	if((((StsVirtualVolumeWizard)wizard).getSeismicVolumes().length == 0)
    		&& (((StsVirtualVolumeWizard)wizard).getVirtualVolumes().length == 0))
    	{
    		mathBtn.setEnabled(false);
    		blendedBtn.setEnabled(false);
    		xplotBtn.setEnabled(false);
    		filterBtn.setEnabled(false);	
    	}
        else
            mathBtn.setSelected(true);
    }
    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout2);

        nameField.initialize(wizard, "volumeName", true, "Name:");
        colorBtn.setEnabled(false);

        typeBox.gbc.anchor = typeBox.gbc.WEST;
        typeBox.add(mathBtn);
        typeBox.add(blendedBtn);
        typeBox.add(xplotBtn);
        typeBox.add(filterBtn);
        typeBox.add(colorBtn);
        typeBox.add(sensorBtn);

        typeGrp.add(mathBtn);
        typeGrp.add(blendedBtn);
        typeGrp.add(xplotBtn);
        typeGrp.add(filterBtn);
        typeGrp.add(colorBtn);
        typeGrp.add(sensorBtn);

        mathBtn.addActionListener(this);
        blendedBtn.addActionListener(this);
        xplotBtn.addActionListener(this);
        filterBtn.addActionListener(this);
        colorBtn.addActionListener(this);
        sensorBtn.addActionListener(this);

        typeDescription.setBackground(Color.lightGray);
        typeDescription.setFont(new java.awt.Font("Dialog", 0, 10));
        typeDescription.setBorder(BorderFactory.createEtchedBorder());
        typeDescription.setText("Creates a virtual volume through the application of math functions between" +
                           " physical volumes, virtual volumes and scalars");
        typeDescription.setLineWrap(true);
        typeDescription.setWrapStyleWord(true);
        panel.add(nameField,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 9, 5, 5), 0, 0));
        panel.add(typeBox,  new GridBagConstraints(0, 0, 1, 2, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.VERTICAL, new Insets(0, 9, 0, 0), 0, 0));
        panel.add(typeDescription,  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 5, 0, 5), 0, 0));
        this.add(panel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == mathBtn)
        {
            typeDescription.setText(
                "Creates a virtual volume through the application of math functions between" +
                " physical volumes, virtual volumes and scalars");
            type = StsVirtualVolume.SEISMIC_MATH;
        }
        if(source == xplotBtn)
        {
            typeDescription.setText(
                "Creates a virtual volume by setting all the samples of a selected volume that" +
                " are included or excluded from a selected crossplot to transparent. Effectively" +
                " the volume is sculpted to the crossplot results producing a volume that shows only" +
                " the samples that match the crossplot results or only the samples that didn't match" +
                " the crossplot results.");
            type = StsVirtualVolume.SEISMIC_XPLOT_MATH;
        }
        if(source == blendedBtn)
        {
            typeDescription.setText(
                "Creates a virtual volume by merging two volumes. The merge is accomplished by " +
                " the application of a logical query against the first volume. Where the query " +
                " results is TRUE the sample remains the same, where FALSE the sample is replaced " +
                " with the value from the second volume.");
            type = StsVirtualVolume.SEISMIC_BLEND;
        }
        if(source == colorBtn)
        {
            typeDescription.setText("Creates a virtual volume by blending three or four volumes into one virtual volume " +
                               " where each volume dictates the red, green, blue and transparency value, or the hue, " +
                               " saturation, intensity and transparency value.");
            type = StsVirtualVolume.RGB_BLEND;
        }
        if(source == filterBtn)
       {
           typeDescription.setText("Creates a virtual volume by applying a user defined smoothing operator to " +
                                   " the selected volume.");
           type = StsVirtualVolume.SEISMIC_FILTER;
       }
        if(source == sensorBtn)
        {
            typeDescription.setText("Creates a virtual volume from sensor data where sample values are " +
                                    " determined by sensor proximity and attributes. Once volume is constructed" + 
                                    " volumetric calculations can be done to determine the total volume of" +
                                    " sensor data. ");
            type = StsVirtualVolume.SENSOR_VOLUME;
        }
    }

    public void setVolumeType(int type)
    {
        this.type = type;
    }
    public int getVolumeType()
    {
        return type;
    }

}
