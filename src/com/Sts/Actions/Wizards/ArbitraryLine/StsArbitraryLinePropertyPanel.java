package com.Sts.Actions.Wizards.ArbitraryLine;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsArbitraryLinePropertyPanel extends StsFieldBeanPanel
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private String name = "";
    private Color stsColor;

    private StsModel model = null;

    StsStringFieldBean lineName = new StsStringFieldBean();
    StsColorComboBoxFieldBean lineColor = new StsColorComboBoxFieldBean();

    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsArbitraryLinePropertyPanel(StsWizard wizard, StsWizardStep wizardStep)
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
        model = wizard.getModel();
        lineName.initialize(this, "name", true, "Name:");
        lineColor.initializeColors(this, "stsColor", "Color:", model.getSpectrum("Basic"));
        lineColor.setValueObject(StsColor.RED);

        this.add(lineName, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(4, 4, 4, 4), 35, 0));
        this.add(lineColor,new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 4, 4, 4), 35, 0));
     }


    public String getName() { return name; }
    public Color getColor() { return stsColor; }
    public void setName(String name)
    {
        this.name = name;
        if(name != "") wizard.enableFinish();
    }
    public void setColor(Color color) { stsColor = color; }
}
