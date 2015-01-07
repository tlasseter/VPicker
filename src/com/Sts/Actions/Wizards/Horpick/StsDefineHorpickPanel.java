package com.Sts.Actions.Wizards.Horpick;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineHorpickPanel extends StsJPanel
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;
    private StsHorpick horpick;

    private String name = "SurfaceName";
    private StsColor stsColor;

    private StsModel model = null;

	StsGroupBox panelBox = new StsGroupBox("Define New Horizon");
    StsStringFieldBean horName = new StsStringFieldBean();
    StsColorComboBoxFieldBean horColorBean = new StsColorComboBoxFieldBean();
    StsButton preferBtn = null;

//    GridBagLayout gridBagLayout1 = new GridBagLayout();
//  JPanel jPanel1 = new JPanel();
//  GridBagLayout gridBagLayout2 = new GridBagLayout();
//  GridBagLayout gridBagLayout3 = new GridBagLayout();

    public StsDefineHorpickPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;

        try
        {
            initialize();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {

//        setLayout(gridBagLayout3);
        model = wizard.getModel();
        horName.initialize(this, "name", true, "Name:");
        StsColor[] colors = model.getSpectrum("Basic").getStsColors();
		stsColor = colors[0];
        horColorBean.initializeColors(this, "stsColor", "Color:", colors);
        preferBtn = new StsButton("Pick Preference", "Set pick preferences", wizard, "pickPreferences");
        gbc.fill = gbc.HORIZONTAL;
		add(panelBox);
        panelBox.gbc.fill = gbc.HORIZONTAL;
		panelBox.add(horName);
		panelBox.add(horColorBean);
		panelBox.gbc.gridwidth = 2;
        panelBox.gbc.fill = gbc.NONE;
        panelBox.gbc.anchor = gbc.SOUTH;
		panelBox.add(preferBtn);

//        jPanel1.setBorder(BorderFactory.createEtchedBorder());
//    jPanel1.setLayout(gridBagLayout2);
//    this.add(jPanel1,   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
//            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 4, 6, 6), 2, 0));
/*
    add(horName,   new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(6, 60, 0, 5), 0, 0));
    add(horColor,    new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 30, 0, 5), 0, 0));
    add(preferBtn,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(8, 62, 3, 55), 35, 0));
*/
     }


    public String getName() { return name; }
    public StsColor getStsColor() { return stsColor; }
    public void setName(String name)
    {
        this.name = name;
        if(name != "") wizard.enableFinish();
    }
    public void setStsColor(StsColor color) { stsColor = color; }
}
