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

public class StsTypeFilterVolumePanel extends JPanel implements ActionListener
{
    StsGroupBox typeBox = new StsGroupBox("Type");
    ButtonGroup typeGrp = new ButtonGroup();

    JRadioButton filterBtn = new JRadioButton("Filter");

    StsStringFieldBean nameField = new StsStringFieldBean();
    JTextArea typeDescription = new JTextArea();
    ButtonGroup typeGroup = new ButtonGroup();

    String name = "FilterName";
    int type = StsVirtualVolume.SEISMIC_FILTER;
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    JPanel panel = new JPanel(new GridBagLayout());

    public StsTypeFilterVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }


    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout2);

        nameField.initialize(this,"name",true,"Name:");

        typeBox.gbc.anchor = typeBox.gbc.WEST;

        typeBox.add(filterBtn);
        typeGrp.add(filterBtn);

        filterBtn.addActionListener(this);

        typeDescription.setBackground(Color.lightGray);
        typeDescription.setFont(new Font("Dialog", 0, 10));
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
        if(source == filterBtn)
        {
            typeDescription.setText("Creates a virtual volume by applying a user defined smoothing operator to " + " the selected volume.");
            type = StsVirtualVolume.SEISMIC_FILTER;
        }
    }

    public void setName(String name)
    {
        this.name = name;
    }
    public void setVolumeType(int type)
    {
        this.type = type;
    }
    public int getVolumeType()
    {
        return type;
    }

    public String getName()
    {
        return name;
    }

}
