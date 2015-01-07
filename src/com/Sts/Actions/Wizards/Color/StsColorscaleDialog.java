package com.Sts.Actions.Wizards.Color;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
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

/** Used for editing a colorscale.  Contains an editable StsColorscalePanel and standard dialog buttons. */

public class StsColorscaleDialog extends JDialog implements WindowListener, ActionListener, ChangeListener
{
    private ButtonGroup btnGroup1 = new ButtonGroup();
    public StsColorscale colorscale = null;
    public float[] dataDist = null;
    public StsModel model = null;
    public JPanel jPanel1 = new JPanel();
    public JPanel buttonPanel = new JPanel();
    public JButton selectBtn = new JButton("Select");
    public JButton newBtn = new JButton("New");
    public JButton exitBtn = new JButton("Ok");
    public JButton cancelBtn = new JButton("Cancel");
	StsColorscalePanel observerColorscalePanel = null;
    StsColorscalePanel colorscalePanel = null;
    GridBagLayout gridBagLayout4 = new GridBagLayout();
    Font defaultFont = new Font("Dialog",0,11);
    boolean success = true;
/*
    public StsColorscaleDialog(String title, StsColorscalePanel observerColorscalePanel)
    {
		super((Frame)null, "Colorscale Editor", false);
		classInitialize(observerColorscalePanel);
		try
		{
			jbInit();
			pack();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
    }
*/
	public StsColorscaleDialog(String title, StsColorscalePanel observerColorscalePanel, Frame frame)
	{
		super(frame, "Colorscale Editor", false);
		initialize(observerColorscalePanel);
		try
		{
			jbInit();
			pack();
            if(frame == null) StsToolkit.centerThisOn(this, observerColorscalePanel);
        }
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
    }

	public StsColorscaleDialog(String title, StsModel model, StsColorscalePanel observerColorscalePanel, Frame frame)
	{
        super(model.win3d, "Colorscale Editor", false);
        // this.setLocationRelativeTo(model.win3d);
        this.model = model;
		initialize(observerColorscalePanel);
		try
		{
			jbInit();
			pack();
            if(frame == null) StsToolkit.centerThisOn(this, observerColorscalePanel);
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

	public StsColorscaleDialog(String title, StsModel model, StsColorscale colorscale)
	{
		super(model.win3d, "Colorscale Editor", false);

		this.colorscale = colorscale;
		this.colorscalePanel = new StsColorscalePanel(colorscale);
		colorscalePanel.setName("Dialog panel");
		try
		{
			jbInit();
			pack();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
    }
	// not needed, histogram data is brought in inside observerColorscalePanel in constructor: get it there TJL 1/31/06
/*
    public void initData(float[] dataDist)
    {
        colorscalePanel.setDataHistogram(dataDist);
    }
*/
    private void jbInit() throws Exception
    {
        this.setModal(false);
        this.setSize(new Dimension(200,500));

        colorscalePanel.setFont(null);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());

        selectBtn.setFont(defaultFont);
        selectBtn.addActionListener(this);
        newBtn.setFont(defaultFont);
        newBtn.addActionListener(this);
        exitBtn.setFont(defaultFont);
        exitBtn.addActionListener(this);
        cancelBtn.setFont(defaultFont);
        cancelBtn.addActionListener(this);

        buttonPanel.setBorder(BorderFactory.createEtchedBorder());
        buttonPanel.setLayout(gridBagLayout4);
        buttonPanel.add(selectBtn,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
                ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        buttonPanel.add(newBtn,  new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0
                ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        buttonPanel.add(cancelBtn,  new GridBagConstraints(2, 0, 1, 1, 1.0, 1.0
                ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        buttonPanel.add(exitBtn,  new GridBagConstraints(3, 0, 1, 1, 1.0, 1.0
                ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));

        jPanel1.setLayout(gridBagLayout4);
        jPanel1.setSize(new Dimension(200,500));
        jPanel1.add(colorscalePanel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
                ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(buttonPanel,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
                ,GridBagConstraints.SOUTH, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));

        this.getContentPane().add(jPanel1, BorderLayout.CENTER);
        this.setResizable(false);
        this.addWindowListener(this);

    }

    private void initialize(StsColorscalePanel observerColorscalePanel)
    {
//        if(model != null && model.win3d != null)
//			model.win3d.getCursor3dPanel().addChangeListener(this);
//		if(observerColorscalePanel != null)
//			colorscale.addChangeListener(this);
		this.observerColorscalePanel = observerColorscalePanel;
		this.colorscale = observerColorscalePanel.getColorscale();
		this.dataDist = observerColorscalePanel.getHistogram();
		observerColorscalePanel.setName("Observer Panel");
        colorscalePanel = new StsColorscalePanel(observerColorscalePanel);
		colorscalePanel.setLabelsOn(true);
		colorscalePanel.setName("Dialog Panel");
    }

    public void actionPerformed(ActionEvent e)
    {
        int i;
        Object source = e.getSource();
        if(source == newBtn)
        {
            StsSpectrumDialog sd = new StsSpectrumDialog("New Spectrum", model, colorscalePanel, true);
            sd.setVisible(true);
            Main.logUsage();
        }
        else if(source == selectBtn)
        {
            StsSpectrumSelect ss = new StsSpectrumSelect(colorscalePanel);
            ss.setVisible(true);
            Main.logUsage();
        }
        else if(source == exitBtn)
        {
            endCurrentAction();
        }
        else if(source == cancelBtn)
        {
            cancelCurrentAction();
        }

    }

    public void endCurrentAction()
    {
        success = true;
		if(model != null)
	    {
			StsActionManager actionManager = model.mainWindowActionManager;
			if (actionManager != null)
				actionManager.endCurrentAction();
			if (model.win3d != null)
				model.win3d.getCursor3dPanel().removeChangeListener(this);
		}
        this.setVisible(false);
    }

    public void cancelCurrentAction()
    {
        success = false;
        model.mainWindowActionManager.endCurrentAction();
        this.setVisible(false);
    }

    public void windowClosing(WindowEvent e)
    {
        endCurrentAction();
    }

    public boolean getSuccess() { return success; }

    public void windowDeactivated(WindowEvent e) { }
    public void windowOpening(WindowEvent e) { }
    public void windowActivated(WindowEvent e) { }
    public void windowDeiconified(WindowEvent e) { }
    public void windowIconified(WindowEvent e) { }
    public void windowClosed(WindowEvent e) { }
    public void windowOpened(WindowEvent e) { }

    public void stateChanged(ChangeEvent e)
    {
        colorscalePanel.repaint();
		observerColorscalePanel.repaint();
        model.win3dDisplayAll();
        Main.logUsageTimer();
    }

//
//  Test for dependent constructor
//
    static public void main(String[] args)
    {
        StsColor[] StsColors =
        {
            new StsColor(0, 0, 0, 255, 255),
            new StsColor(255, 0, 255, 255, 255)
        };

        try
        {
            StsModel model = new StsModel();
            StsSpectrum spectrum = StsSpectrum.constructor("Default Spectrum", StsColors, 255);
            StsColorscale colorscale = new StsColorscale("Test", spectrum, 0.0f, 255.0f);
            StsColorscaleDialog dialog = new StsColorscaleDialog("Test", model, colorscale);
            dialog.setVisible(true);
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
}
