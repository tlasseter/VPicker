package com.Sts.Workflow;

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsWorkflowSelectPanel extends JPanel implements ActionListener
{
    StsWin3d win3d = null;

    ImageIcon vendorIcon = null;
    ImageIcon imageIcon = null;
    Image image;
    
    JButton s2sLogoBtn = new JButton();
    JButton vendorLogoBtn = new JButton();

    JLabel label1, label2;
    
    public StsWorkflowSelectPanel(StsWin3d win3d)
    {
        try
        {
            this.win3d = win3d;
            if(win3d.model.workflowPlugInNames == null) return;
            if(win3d.model.workflowPlugInNames.length == 1)
                Main.singleWorkflow = true;
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    void jbInit() throws Exception
    {
    	this.setBackground(Color.BLACK);
        image = StsIcon.createImage(Main.vendorName + "SplashScreen.jpg");
    }

    public void paintComponent(Graphics g)
    {
    	super.paintComponent(g); 
        if(image != null) g.drawImage(image, 0,0,this.getWidth(),this.getHeight(),this);
    }
    
    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == s2sLogoBtn)
        {
            if(!StsToolkit.launchWebPage("http://www." + Main.vendorName + ".com"))
                new StsMessage(win3d,StsMessage.WARNING,"Unable to find or launch webpage: http://www." + Main.vendorName + ".com,\n"
                              + " possibly disconnected from Internet or a browser is not in your system path");
        }
        else if(source == vendorLogoBtn)
        {
            if(!StsToolkit.launchWebPage(vendorLogoBtn.getToolTipText()))
                new StsMessage(win3d,StsMessage.WARNING,"Unable to find webpage: " + vendorLogoBtn.getToolTipText() + ", possibly disconnected from Internet");
        }
    }

}
