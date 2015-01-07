package com.Sts.Actions.Wizards.WizardHeaders;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.jnlp.*;
import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsHeaderPanel extends JPanel implements ActionListener
{
    private JLabel headerTxt = new JLabel();
    private ButtonGroup zGroup = new ButtonGroup();
    private JLabel subtitleTxt = new JLabel();
    private JTextArea infoTxt = new JTextArea();
    private JButton iconBtn = new JButton();
    private String logoLink = "http://www.s2ssystems.com";
    private TitledBorder titledBorder1;
    private JPanel jPanel1 = new JPanel();
    private JPanel jPanel2 = new JPanel();
    private GridBagLayout gridBagLayout2 = new GridBagLayout();
    private GridBagLayout gridBagLayout1 = new GridBagLayout();

    private DownloadService downloadService;
    GridBagLayout gridBagLayout3 = new GridBagLayout();
    JButton infoBtn = new JButton();
    boolean infoOn = false;
    private StsWizardDialog dialog = null;

    public static int displayableHeight = 72;
    public StsHeaderPanel()
    {
        try
        {
            jbInit();
            setDefaultLogoAndLink();
            hideInfoText();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }
    public void setDialog(StsWizardDialog dialog)
    {
        this.dialog = dialog;
    }

    void jbInit() throws Exception
    {
        titledBorder1 = new TitledBorder("");
        headerTxt.setFont(new java.awt.Font("Dialog", 1, 14));
        headerTxt.setHorizontalAlignment(SwingConstants.CENTER);
        headerTxt.setText("Title");
        subtitleTxt.setFont(new java.awt.Font("Dialog", 0, 11));
        subtitleTxt.setHorizontalAlignment(SwingConstants.CENTER);
        subtitleTxt.setHorizontalTextPosition(SwingConstants.CENTER);
        subtitleTxt.setText("subTitle");
        iconBtn.setAlignmentY((float) 0.0);
        iconBtn.setBorder(null);
        iconBtn.setMaximumSize(new Dimension(72, displayableHeight));
        iconBtn.setPreferredSize(new Dimension(72, displayableHeight));
        iconBtn.setBorderPainted(false);
        iconBtn.setHorizontalTextPosition(SwingConstants.CENTER);
        iconBtn.setMargin(new Insets(0, 0, 0, 0));
        iconBtn.setMnemonic('0');
        iconBtn.setToolTipText(logoLink);
        iconBtn.addActionListener(this);

        infoTxt.setFont(new java.awt.Font("Dialog", 0, 10));
        infoTxt.setLineWrap(true);
        infoTxt.setEditable(false);
        infoTxt.setWrapStyleWord(true);
        infoTxt.setBackground(SystemColor.menu);
        infoTxt.setRequestFocusEnabled(true);
        infoTxt.setText("None Available");
        jPanel2.setBackground(SystemColor.menu);
        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setLayout(gridBagLayout3);

        this.setLayout(gridBagLayout2);
        jPanel1.setLayout(gridBagLayout1);
        //infoBtn.setBackground(SystemColor.black);
        infoBtn.setAlignmentY((float) 0.0);
        infoBtn.setMaximumSize(new Dimension(15, displayableHeight));
        infoBtn.setPreferredSize(new Dimension(15, displayableHeight));
        infoBtn.setToolTipText("Show instructions on how to use this screen");
        infoBtn.setIcon(createIcon("Info.gif"));
        infoBtn.setMargin(new Insets(0, 0, 0, 0));
        infoBtn.addActionListener(this);
        
        jPanel1.add(subtitleTxt,             new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 10, 0), 10, 0));
        jPanel1.add(headerTxt,         new GridBagConstraints(2, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 0), 0, 10));
        jPanel1.add(iconBtn,  new GridBagConstraints(1, 0, 1, 2, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 3, 0, 0), 0, 0));
        jPanel1.add(infoBtn,  new GridBagConstraints(0, 0, 1, 2, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 3, 0, 0), 0, 0));

        this.add(jPanel1,      new GridBagConstraints(0, 0, 1, 1, 0.5, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 2, 4), 0, 0));

        jPanel2.add(infoTxt,              new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 5, 5), 0, 0));
        this.add(jPanel2,       new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 0, 5), 0, 0));
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if((source == iconBtn) && (logoLink != null))
        {
            if(!StsToolkit.launchWebPage(logoLink))
            {
                StsMessageFiles.infoMessage("Unable to load webpage: " + logoLink + ", trying default page.");
                if(!StsToolkit.launchWebPage("http://www.s2ssystems.com"))
                {
                    StsMessageFiles.infoMessage("Unable to access Internet, failed to load default page.");
                    new StsMessage(new Frame(), StsMessage.WARNING, logoLink + " not accessible, possibly disconnected from Internet");
                }
            }
        }
        else if(source == infoBtn)
        {
            toggleInfo();
        }
    }

    public void setTitle(String title)  { headerTxt.setText(title); }
    public void setSubtitle(String subtitle)  { subtitleTxt.setText(subtitle); }
    public void addInfoLine(StsWizardDialog dialog, String line)
    {
        this.dialog = dialog;
        String info = infoTxt.getText();
        if(info.equals("None Available"))
            info = line;
        else
            info = info + "\n" + line;
        infoTxt.setText(info);
    }
    public void setInfoText(StsWizardDialog dialog, String info)
    {
        this.dialog = dialog;
        infoTxt.setText(info);
    }
    public void hideInfoText()
    {
        this.remove(jPanel2);
        this.setSize(this.getWidth(), this.getHeight()-jPanel2.getHeight());
        repaint();      
        if(dialog != null)
            dialog.pack();

        infoOn = false;
    }
    public void showInfoText()
    {     
        this.add(jPanel2,       new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 0, 5), 0, 0));
        this.setSize(this.getWidth(), this.getHeight()+jPanel2.getHeight());

        repaint();
        if(dialog != null)
            dialog.pack();

        infoOn = true;
    }
    private void toggleInfo()
    {
        if(infoOn)
            hideInfoText();
        else
            showInfoText();
        
        resizeDialog();
    }

    private void resizeDialog()
    {
        Dimension dim = dialog.getSize();    	
    	if(infoOn)
    	{
            dim.height = dim.height + jPanel2.getHeight();
    	}
    	else
    	{
    		dim.height = dim.height - jPanel2.getHeight();
    	}
        dialog.setPreferredSize(dim);    	
        repaint();
        if(dialog != null)
            dialog.pack();    	
    }
    
    public void setDefaultLogoAndLink()
    {
        iconBtn.setIcon(createIcon(Main.vendorName + "Logo.gif"));
        this.logoLink = "www." + Main.vendorName + ".com";
        iconBtn.setToolTipText(logoLink);
    }

    public void setLogo(String imageName)
    {
        iconBtn.setIcon(createIcon(imageName));
    }

    public void setLink(String logoLink)
    {
        this.logoLink = logoLink;
        iconBtn.setToolTipText(logoLink);
    }

    protected ImageIcon createIcon(String name)
    {
        return StsIcon.createIcon(getClass(), name);
    }

}
