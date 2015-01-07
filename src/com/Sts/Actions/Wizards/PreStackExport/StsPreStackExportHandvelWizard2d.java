package com.Sts.Actions.Wizards.PreStackExport;

import com.Sts.MVC.*;

public class StsPreStackExportHandvelWizard2d extends StsPreStackExportHandvelWizard 
{

	public StsPreStackExportHandvelWizard2d(StsActionManager actionManager)
	{
		super(actionManager);
	}

	@Override
	protected StsExportHandvel[] createHandVelSteps() 
	{
		return new StsExportHandvel[] { new StsExportHandvel2d(this) };
	}

}
