package com.Sts.UI;

import com.Sts.Interfaces.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

public class StsListScrollPanel extends JPanel implements Runnable
{
	private int type;
	private ListPanel listPanel;
	private JScrollBar vsb = new JScrollBar(JScrollBar.VERTICAL);
	private JScrollBar hsb = new JScrollBar(JScrollBar.HORIZONTAL);
	private Border border = BorderFactory.createEtchedBorder();

	public StsListScrollPanel(StsListModel listModel, int type)
	{
		this.type = type;
		vsb.setBorder(border);
		hsb.setBorder(border);
		setLayout(new BorderLayout());
		listPanel = new ListPanel(listModel);
		listPanel.setBackground(Color.white);
		add(vsb, BorderLayout.EAST);
		add(hsb, BorderLayout.SOUTH);
		add(listPanel, BorderLayout.CENTER);

		vsb.addAdjustmentListener(new AdjustmentListener()
		{
			public void adjustmentValueChanged(AdjustmentEvent e)
			{
				JScrollBar sb = (JScrollBar)e.getSource();
				listPanel.setTopIndexByPixelValue(e.getValue());
				repaint();
			}
		});
		hsb.addAdjustmentListener(new AdjustmentListener()
		{
			public void adjustmentValueChanged(AdjustmentEvent e)
			{
				JScrollBar sb = (JScrollBar)e.getSource();
				listPanel.setLeftPixelValue(e.getValue());
				adjustHorizontalScroll();
				repaint();
			}
		});

		final Runnable r = this;
		listModel.addListDataListener(new ListDataListener()
		{

			public void intervalAdded(ListDataEvent e)
			{
                StsToolkit.runLaterOnEventThread(r);
			}

			public void intervalRemoved(ListDataEvent e)
			{
			}

			public void contentsChanged(ListDataEvent e)
			{
			}
		});

		addComponentListener(new ComponentAdapter()
		{
			public void componentResized(ComponentEvent e)
			{
				adjustVerticalScroll();
				adjustHorizontalScroll();
			}
		});
	}

	public void run()
	{
		adjustVerticalScrollToBottom();
		repaint();
	}

	private void adjustVerticalScroll()
	{
		int max = listPanel.getPreferredSize().height;
		int extent = listPanel.getSize().height;
		int value = Math.max(0, Math.min(hsb.getValue(), max - extent));
		adjustVerticalScroll(value, extent, max);
	}

	private void adjustVerticalScroll(int value, int extent, int max)
	{
		listPanel.setTopIndexByPixelValue(Math.max(0, max - extent));
		vsb.setVisible(extent < max);
		vsb.setUnitIncrement(listPanel.getUnitHeight());
		vsb.setBlockIncrement(extent - vsb.getUnitIncrement());
		vsb.setValues(value, extent, 0, max);
	}

	private void adjustVerticalScrollToBottom()
	{
		int max = listPanel.getPreferredSize().height;
		int extent = listPanel.getSize().height;
		int value = Math.max(0, max - extent);
		adjustVerticalScroll(value, extent, max);
	}

	private void adjustHorizontalScroll()
	{
		int max = listPanel.getPreferredSize().width;
		int extent = listPanel.getSize().width;
		int value = Math.max(0, Math.min(hsb.getValue(), max - extent));
		hsb.setVisible(extent < max);
		hsb.setUnitIncrement(listPanel.getUnitWidth());
		hsb.setBlockIncrement(Math.max(0, Math.min((max - extent) / 4, hsb.getUnitIncrement())));
		hsb.setValues(value, extent, 0, max);
	}

	public void paint(Graphics g)
	{
		super.paint(g);
	}

}

class ListPanel extends JPanel
{
	private StsListModel listModel;
	private int topIndex = 0;
	private int leftPixelValue = 0;
	private FontMetrics fm = null;
	private float floatFontWidth;
	private int fontHeight, fontWidth;

	// used to get average font width
	static final String widthString = "abcdefghijklmnopqrstuvwxyz1234567890.ABCDEFG";

	public ListPanel()
	{
	}

	public ListPanel(StsListModel listModel)
	{
		this.listModel = listModel;
	}

	public void setListModel(StsListModel listModel)
	{
		this.listModel = listModel;
	}

	public void paintComponent(Graphics g)
	{
		Color color = g.getColor();
		super.paintComponent(g);
		g.setColor(color);

		Dimension size = getSize();
		Insets insets = getInsets();
		int y = insets.top + fontHeight;

		if(listModel == null)return;

		int length = listModel.getSize();
		for(int i = topIndex; i < length; ++i, y += fontHeight)
		{
			String string = (String)listModel.getElementAt(i);
			if(string != null)g.drawString(string, leftPixelValue + fontWidth / 2, y);
			if(y + fontHeight > size.height - insets.bottom)break;
		}
	}

	public void setTopIndexByPixelValue(int pixelValue)
	{
		topIndex = pixelValue / fontHeight;
	}

	public void setLeftPixelValue(int pixelValue)
	{
		leftPixelValue = -pixelValue;
	}

	public int getUnitHeight()
	{
		return fontHeight;
	}

	public int getUnitWidth()
	{
		return 5 * fontWidth;
	}

	public Dimension getPreferredSize()
	{
		if(listModel == null)return new Dimension(100, 10);

		if(fm == null)computeFontDimensions();
		Dimension dim = new Dimension();
		int listLength = listModel.getSize();
		dim.height = fontHeight * (listLength + 1);
		String longestString = listModel.getLongestString();
		dim.width = fm.stringWidth(longestString) + fontWidth;
		return dim;
	}

	private void computeFontDimensions()
	{
		Graphics g = getGraphics();

		try
		{
			fm = g.getFontMetrics();
			fontHeight = fm.getHeight();
//			int listLength = listModel.getSize();
			int stringLength = widthString.length();
			int width = fm.stringWidth(widthString);
			floatFontWidth = (float)width / stringLength;
			fontWidth = (int)floatFontWidth;
			if(floatFontWidth - fontWidth > 0.5f)fontWidth++;
		}
		finally
		{
			g.dispose();
		}
	}
}
