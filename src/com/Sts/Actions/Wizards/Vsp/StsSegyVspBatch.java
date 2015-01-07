package com.Sts.Actions.Wizards.Vsp;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsSegyVspBatch extends StsPostStackBatch
{
	public StsSegyVspBatch(StsWizard wizard)
	{
		super(wizard);
        header.setTitle("VSP Volumes Set to Process");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Vsp");                
    }
}
