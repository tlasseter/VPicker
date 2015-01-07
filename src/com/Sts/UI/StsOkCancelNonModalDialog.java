package com.Sts.UI;

import com.Sts.Interfaces.*;

import java.awt.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author TJLasseter
 * @version c63
 */
abstract public class StsOkCancelNonModalDialog extends StsOkCancelDialog
{
	public StsOkCancelNonModalDialog(Frame frame, StsDialogFace okCancelObject, String title)
	{
		super(frame, okCancelObject, title, false);
	}

	abstract public StsOkCancelNonModalDialog getDialog();
	abstract public void deleteDialog();

	public void close()
	{
		super.close();
		deleteDialog();
	}

	public void ok()
	{
		super.ok();
		deleteDialog();
	}

	public void cancel()
	{
		super.cancel();
		deleteDialog();
	}
}
