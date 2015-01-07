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

public class StsSelectHorpickPanel extends StsJPanel
{
    private StsHorpickWizard wizard;
    private StsSelectHorpick wizardStep;

    private StsModel model = null;
    private StsHorpick selectedHorpick = null;

	StsGroupBox groupBox = null;
	StsListFieldBean horpickList = new StsListFieldBean();
	StsButton newHorpickButton = new StsButton();

    public StsSelectHorpickPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsHorpickWizard)wizard;
        this.wizardStep = (StsSelectHorpick)wizardStep;
    }

    public void initialize()
    {
        if(groupBox != null)
        {
            remove(groupBox);
            constructPanel();
            wizard.rebuild();
        }
        else
            constructPanel();
    }

    private void constructPanel()
    {
        groupBox = new StsGroupBox("Select or create pick surface");
        model = wizard.getModel();
		StsHorpickWizard horpickWizard = wizard;
        StsObject[] horpicks = model.getObjectList(StsHorpick.class);
        groupBox.gbc.fill = groupBox.gbc.BOTH;
		if(horpicks.length > 0)
		{
			horpickList.initialize(horpickWizard, "selectedHorpick", "Existing Horiz Pick", horpicks);
			groupBox.add(horpickList);
		}
        groupBox.gbc.fill = groupBox.gbc.NONE;
        groupBox.gbc.anchor = groupBox.gbc.SOUTH;
		groupBox.gbc.gridwidth = 2;
		newHorpickButton.initialize("Pick New Surface", "Interactively pick a new surface.", horpickWizard, "createNewHorpick", null);
		groupBox.add(newHorpickButton);
        gbc.fill = gbc.HORIZONTAL;
		add(groupBox);
    }
}
