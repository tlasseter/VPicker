package com.Sts.Actions.Wizards.Velocity;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: jS2S development</p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author T.J.Lasseter
 * @version c51c
 */
public class StsFromSonicsPanel extends StsJPanel
{
	StsModel model;
	StsBooleanFieldBean fromSonicsBean;
	StsBooleanFieldBean fromVelfBean;
	boolean useSonics = true;
	boolean useVelf = false;
    private StsGroupBox importGroupBox = new StsGroupBox("Import Profiles");

    private StsGroupBox defineGroupBox = new StsGroupBox("Define Volume");

    StsSeisVelWizard wizard = null;
    StsWizardStep wizardStep;

    static final boolean debug = false;

    public StsFromSonicsPanel(StsFromSonics assignVelocity)
    {
		super();
        wizard = (StsSeisVelWizard) assignVelocity.getWizard();
        wizardStep = assignVelocity;
		model = assignVelocity.getModel();
        constructPanel();
    }

    private void constructPanel()
    {
        importGroupBox.gbc.fill = gbc.HORIZONTAL;
        //importGroupBox.addEndRow(surfaceListBean);
        //gbc.fill = gbc.HORIZONTAL;
        //gbc.gridwidth = 2;
        //this.addEndRow(importGroupBox);
		fromSonicsBean = new StsBooleanFieldBean(this, "fromSonics", true, "From Well T/D functions", true);
		fromVelfBean = new StsBooleanFieldBean(this, "fromVelf", false, "From Additional Imported T/D functions", true);


		 gbc.gridwidth = 2;
		 gbc.fill = gbc.HORIZONTAL;

        defineGroupBox.gbc.fill = gbc.HORIZONTAL;
        //defineGroupBox.addEndRow(zIncrementBean);
		defineGroupBox.add(fromSonicsBean);
		defineGroupBox.add(fromVelfBean);
        this.addEndRow(defineGroupBox);
    }

	public void setFromSonics(boolean value)
		{
			useSonics = value;
		}
    public boolean getFromSonics() { return useSonics; }
	public void setFromVelf(boolean value)
	{
		useVelf = value;
	}
    public boolean getFromVelf() { return useVelf; }


}