package com.Sts.UI.Beans;

/**
 * <p>Title: Field Beans Development</p>
 * <p>Description: General beans for generic panels.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.0
 */

import java.awt.*;

public class StsGroupBox extends StsJPanel
{
	protected String label;

	public StsGroupBox()
	{
		this("");
	}

	public StsGroupBox(String label)
	{
		this(label, true);
	}

	public StsGroupBox(String label, boolean addInsets)
	{
		initialize(label, addInsets);
	}

	static public StsGroupBox noInsets(String label)
	{
		return new StsGroupBox(label, false);
	}

    public void initialize(String label, boolean addInsets)
    {
        this.label = label;
        initializeLayout(addInsets);
    }

	public Dimension getMinimumSize()
	{
		return getPreferredSize();
	}

	public void setLabel(String label)
	{
		this.label = label;
		repaint();
	}

	/** @DEPRECATED */
	public Insets insets()
	{
		Font f = getFont();
		int inset = 0;
		if(f != null)
			inset = getFontMetrics(f).getHeight();
		return new Insets(inset, inset, inset, inset);
	}

	public void paintComponent(Graphics g)
	{
		super.paintComponent(g);
		// Save the old font
		Font of = g.getFont();
		Font f = getFont();
		g.setFont(f);
		FontMetrics fm = g.getFontMetrics(f);
		int ascent = fm.getAscent();

		int ww = ascent / 2 + 1;

		Rectangle r = getBounds();
		if(r == null)
			return;
		int xs = ww;
		int ys = ww;
		int xe = r.width - ww - 1;
		int ye = r.height - ww - 1;

		g.setColor(getForeground());
		g.drawString(label, xs + xs, ascent);
		int tw = fm.stringWidth(label);

		for(int i = 1; i >= 0; i--)
		{
			g.setColor(i == 0 ? SystemColor.controlShadow : SystemColor.controlLtHighlight);
			g.drawLine(xs + xs - 1, ys + i, xs + i, ys + i); // top left piece
			g.drawLine(xs + i, ys + i, xs + i, ye + i); // left
			g.drawLine(xs + i, ye + i, xe + i, ye + i); // bottom
			g.drawLine(xe + i, ye + i, xe + i, ys + i); // right
			g.drawLine(xe + i, ys + i, xs + xs + tw, ys + i); // top right piece
		}
		// Restore the old font
		g.setFont(of);
	}

	public static void main(String[] args)
	{
	}
}
