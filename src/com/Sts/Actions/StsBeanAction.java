package com.Sts.Actions;

import com.Sts.MVC.StsActionManager;
import com.Sts.UI.Beans.StsFieldBean;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: 5/23/11
 * Time: 7:47 AM
 * To change this template use File | Settings | File Templates.
 */
abstract public class StsBeanAction extends StsAction
{
	protected StsFieldBean actionBean;

	public StsBeanAction() { }

	public StsBeanAction(StsActionManager actionManager, StsFieldBean actionBean)
	{
		super(actionManager);
		this.actionBean = actionBean;
	}

	public StsBeanAction(StsActionManager actionManager, StsFieldBean actionBean, boolean canInterrupt)
	{
		super(actionManager, canInterrupt);
		this.actionBean = actionBean;
	}
}
