package com.Sts.Actions.Wizards.PreStackExport;

import com.Sts.MVC.*;

public class StsPreStackExportHandvelWizard3d extends StsPreStackExportHandvelWizard 
{

	public StsPreStackExportHandvelWizard3d(StsActionManager actionManager)
	{
		super(actionManager);
	}

	@Override
	protected StsExportHandvel[] createHandVelSteps() 
	{
		return new StsExportHandvel[] { new StsExportHandvel3d(this) };
	}

}
