package com.Sts.Actions.Wizards.LogEdit;

import com.Sts.DBTypes.StsFracture;
import com.Sts.DBTypes.StsFractureSet;
import com.Sts.DBTypes.StsObjectRefList;
import com.Sts.DBTypes.StsSectionEdge;
import com.Sts.IO.StsAsciiFile;
import com.Sts.IO.StsFile;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.StsProject;
import com.Sts.Types.StsPoint;
import com.Sts.Types.StsRotatedGridBoundingBox;
import com.Sts.UI.Beans.StsBooleanFieldBean;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Beans.StsStringFieldBean;
import com.Sts.UI.StsButton;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsYesNoDialog;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsParameters;
import com.Sts.Utilities.StsStringUtils;

import javax.swing.*;
import java.awt.*;
import java.io.File;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsNewLogDialog extends JDialog
{
    StsModel model;

	StsJPanel panel = StsJPanel.addInsets();
	StsGroupBox selectBox;
	StsStringFieldBean nameBean;

	StsJPanel buttonPanel = StsJPanel.addInsets();
	StsButton okButton = new StsButton("Process", "Create new log curve.", this, "ok");
	StsButton cancelButton = new StsButton("Cancel", "Cancel this operation.", this, "cancel");

    public final static byte OK = 0;
    public final static byte CANCELED = 1;
    byte mode = OK;

    protected StsNewLogDialog(StsModel model, Frame frame, String title, boolean modal)
    {
        super(frame, title, modal);
        this.model = model;
        constructBeans();
        constructPanel();
    }

    protected void constructBeans()
    {
        selectBox = new StsGroupBox("Export fracture set");
        nameBean = new StsStringFieldBean(this, "logType", "New log type name");
    }

    protected void constructPanel()
    {

	    this.getContentPane().add(panel);
		this.setTitle("Well Export Parameters");
		panel.add(selectBox);
		panel.add(buttonPanel);
		selectBox.add(nameBean);

	    buttonPanel.addToRow(okButton);
		buttonPanel.addEndRow(cancelButton);
    }

    public void ok()
	{
		mode = OK;
        setVisible(false);
	}

	public void cancel()
	{
		mode = CANCELED;
		setVisible(false);
 	}

    public byte getMode() { return mode; }

    static public void main(String[] args)
    {

    }
}
