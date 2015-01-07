package com.Sts.Actions.Wizards.PostStack3d;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsPostStackTraceDefinitionPanel3d extends StsPostStackTraceDefinitionPanel
{
    
    public StsPostStackTraceDefinitionPanel3d(StsSeismicWizard wizard, StsPostStackTraceDefinition wizardStep)
    {
        super( wizard, wizardStep);
        this.wizardStep = wizardStep;
    }

    public void buildTablePanel()
    {
        String[] columnNames = {"stemname", "xOrigin", "yOrigin", "rowNumMin", "rowNumMax", "colNumMin", "colNumMax", "statusString"};
        String[] columnTitles = {"Name", "Start X", "Start Y", "Start Inline", "End Inline", "Start XLine", "End XLine", "Status"};
        StsEditableTableModel tableModel = new StsEditableTableModel(StsSeismicBoundingBox.class, columnNames, columnTitles, true);
        volumeStatusTablePanel = new StsTablePanelNew(tableModel);
        volumeStatusTablePanel.setLabel("Volumes");
        volumeStatusTablePanel.setSize(400, 100);
        volumeStatusTablePanel.initialize();
    }

    static public void main(String[] args)
	{
		JFrame frame = new JFrame("Test Panel");
		frame.setSize(300, 200);

		Container contentPane = frame.getContentPane();
		StsPostStackTraceDefinitionPanel3d panel = new StsPostStackTraceDefinitionPanel3d(null, null);
		contentPane.add(panel);
		StsToolkit.centerComponentOnScreen(frame);
		frame.pack();
		frame.setVisible(true);
	}
}
