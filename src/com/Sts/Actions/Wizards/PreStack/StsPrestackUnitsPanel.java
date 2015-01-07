package com.Sts.Actions.Wizards.PreStack;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.PreStack3d.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */
public class StsPrestackUnitsPanel extends StsPoststackUnitsPanel
{
	private StsBooleanFieldBean nmoChk;

	public StsPrestackUnitsPanel(StsSeismicWizard wizard)
	{
        super(wizard);
    }

	public void initialize()
	{
		super.initialize();
		nmoChk.getValueFromPanelObject();
	}

	protected void constructBeans()
	{
        super.constructBeans();
		nmoChk = new StsBooleanFieldBean(this, "NMOed", "NMOed", false);
        nmoChk.setToolTipText("Is the data NMOed?");
	}

	protected void constructPanel()
	{
        super.constructPanel();
        gbc.gridwidth = 2;
		addEndRow(nmoChk);
	}

    public boolean getNMOed()
    {
        return ((StsPreStackWizard)wizard).getNMOed();
    }

    public void setNMOed(boolean value)
    {
        ((StsPreStackWizard)wizard).setNMOed(value);
    }
}
