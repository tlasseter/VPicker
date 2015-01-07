package com.Sts.UI;

import com.Sts.MVC.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.Progress.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsSplashScreen extends JDialog
{
	private JPanel panel = new JPanel();
	private JLabel s2sIcon = new JLabel();
	public StsProgressBar progressBar = StsProgressBar.constructor();
	private static StsSplashScreen splash = null;
	static private String version = "Unknown";
	private static String vendor = "GeoCloudRealTime";

    static public StsSplashScreen instance = null;

    private StsSplashScreen(String version, String vendor)
	{
		super();
		this.version = version;
		this.vendor = vendor;
		try
		{
			jbInit();
			this.setLocationRelativeTo(null);
			setModal(false);
			setVisible(true);
			toFront();
			repaint();
            instance = this;
        }
		catch (Exception exception)
		{
			exception.printStackTrace();
		}
	}

	public static StsSplashScreen createAndShow(String version)
	{
		if (SwingUtilities.isEventDispatchThread())
		{
			splash = new StsSplashScreen(version, Main.vendorName);
			splash.show();
		}
		else
		{
			Runnable showSplash = new Runnable()
			{
				public void run()
				{
					splash = new StsSplashScreen("", Main.vendorName);
                    splash.paintImmediately();
                }
			};
			try
			{
				SwingUtilities.invokeAndWait(showSplash);
			}
			catch (Exception ex)
			{
				ex.printStackTrace();
			}
		}
		return splash;
	}

	static public String getVendor()
	{
		return vendor;
	}

	private void jbInit() throws Exception
	{
		this.setUndecorated(true);
		//this.setAlwaysOnTop(true);
		Icon icon = StsIcon.createIcon(vendor + "SplashLogo.gif");
		if (icon == null)
			icon = StsIcon.createIcon("s2ssystemsSplashLogo.gif");
		int height = icon.getIconHeight() + 76;
		this.setSize(320, height);

		panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
		//jPanel1.setLayout(borderLayout2);
		panel.setBackground(Color.white);
		panel.setBorder(BorderFactory.createLineBorder(Color.black));
		s2sIcon.setBackground(new Color(252, 250, 252));
		s2sIcon.setBorder(BorderFactory.createEmptyBorder(12, 0, 4, 0));
		s2sIcon.setHorizontalAlignment(SwingConstants.CENTER);
		s2sIcon.setHorizontalTextPosition(SwingConstants.CENTER);
		s2sIcon.setIcon(icon);
		s2sIcon.setText("Copyright(c) S2S Systems LLC 2001-2011");
		s2sIcon.setVerticalTextPosition(SwingConstants.BOTTOM);
		JLabel versionLbl = new JLabel("Version: " + StsModel.version);
		s2sIcon.setAlignmentX(Component.CENTER_ALIGNMENT);
		versionLbl.setAlignmentX(Component.CENTER_ALIGNMENT);
		panel.add(s2sIcon, BorderLayout.CENTER);
		panel.add(versionLbl);
		this.getContentPane().add(panel, java.awt.BorderLayout.CENTER);
		progressBar.setAlignmentX(Component.CENTER_ALIGNMENT);
		panel.add(progressBar);
	}

	public void setValue(double value)
	{
		progressBar.setValue((int)(value * 100));
	}

    public void setStringImmediate(String progressDescription)
	{
        progressBar.setStringImmediate(progressDescription);
	}

    public void initialize(int maxValue)
    {
        progressBar.initialize(maxValue);
    }

    public void initializeImmediate(int maxValue)
    {
        progressBar.initializeImmediate(maxValue);
    }

    public void setValueAndStringImmediate(int value, String string)
    {
        progressBar.setValueImmediate(value);
        progressBar.setStringImmediate(string);
        paintImmediately();
    }

    public void setValueImmediate(int value)
    {
        progressBar.setValueImmediate(value);
        paintImmediately();
    }

    public void paintImmediately()
	{
        panel.paintImmediately(0, 0, panel.getWidth(), panel.getHeight());
	}

    public void setDescriptionAndLevel(String progressDescription, int level)
	{
        progressBar.setDescriptionAndLevel(progressDescription, level);
	}

    public void finished()
    {
        progressBar.finished();
    }
}
