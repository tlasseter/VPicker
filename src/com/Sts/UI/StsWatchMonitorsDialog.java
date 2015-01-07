package com.Sts.UI;

import com.Sts.MVC.*;

import javax.swing.*;

public class StsWatchMonitorsDialog extends JDialog
{
    StsModel model;
    StsWatchMonitorPanel panel;
    
	public StsWatchMonitorsDialog(StsModel model, boolean modal, boolean hasNewBtn)
	{
		super();
		setTitle("Watch Monitors");
        this.model = model;
        panel = new StsWatchMonitorPanel(model, this, hasNewBtn);
        panel.initialize();
        this.add(panel);
        pack();
        validate();
	}	
	
	public void updateDialog()
	{
		panel.updateTable();
	}
}
