package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.Utilities.*;

public class StsSpectrumClass extends StsClass implements StsSerializable
{
    static StsColor[] defaultSeismicSpectrumColors = new StsColor[]
    {
        new StsColor(0, 128, 0, 0, 255),
        new StsColor(17, 255, 64, 0, 255),
        new StsColor(35, 255, 191, 0, 255),
        new StsColor(52, 255, 255, 64, 255),
        new StsColor(70, 255, 255, 128, 255),
        new StsColor(87, 255, 255, 191, 255),
        new StsColor(105, 246, 246, 246, 255),
        new StsColor(122, 255, 255, 255, 255),
        new StsColor(140, 228, 228, 228, 255),
        new StsColor(157, 191, 191, 191, 255),
        new StsColor(175, 155, 155, 155, 255),
        new StsColor(193, 118, 118, 118, 255),
        new StsColor(210, 82, 82, 82, 255),
        new StsColor(228, 46, 46, 46, 255),
        new StsColor(254, 9, 9, 9, 255)
    };

    static StsColor[] defaultSimulationSpectrumColors = new StsColor[]
    {
        new StsColor(0, 255, 0, 0, 255),
        new StsColor(128, 0, 255, 0, 255),
        new StsColor(254, 0, 0, 255, 255)
    };

    static StsColor[] defaultPropertySpectrumColors = new StsColor[]
    {
        new StsColor(0, 255, 0, 0, 255),
        new StsColor(128, 0, 0, 0, 255),
        new StsColor(254, 0, 0, 255, 255)
    };

    static StsColor[] defaultTempColors = new StsColor[]
    {
       new StsColor(0,  0,  0,    255, 255),
       new StsColor(128,  0,   255,  0,   255),
       new StsColor(254, 255, 0,  0,   255)
    };

    static StsColor[] defaultSemblanceColors = new StsColor[]
    {
    	new StsColor(0,   0,     0,    255, 255),
    	new StsColor(57,  0,   255,    255, 255),
    	new StsColor(99,  0,   255,      0, 255),
    	new StsColor(142, 255, 255,      0, 255),
    	new StsColor(184, 255,   0,      0, 255),
    	new StsColor(254, 255,   0,    255, 255)

//       new StsColor(0,   0,   0,    255, 255),
//       new StsColor(84,  0,   255,  0,   255),
//       new StsColor(168, 255, 255,  0,   255),
//       new StsColor(254, 180, 0,    0,   255)
       
 /*
       new StsColor(0 , 180, 0, 0, 0),
       new StsColor(84 ,255, 255, 0, 255),
       new StsColor(168, 0, 255, 0, 255),
       new StsColor(254, 0, 0, 255, 255)
 */
    };

    static StsColor[] defaultWSatColors = new StsColor[]
    {
    	new StsColor(0,   255,   0,    255, 255),
    	new StsColor(57,  255,   0,      0, 255),
    	new StsColor(99,  255, 255,      0, 255),
    	new StsColor(142, 0,   255,      0, 255),
    	new StsColor(184, 0,   255,    255, 255),
    	new StsColor(254, 0,     0,    255, 255)

//       new StsColor(0,   0,   0,    255, 255),
//       new StsColor(84,  0,   255,  0,   255),
//       new StsColor(168, 255, 255,  0,   255),
//       new StsColor(254, 180, 0,    0,   255)

 /*
       new StsColor(0 , 180, 0, 0, 0),
       new StsColor(84 ,255, 255, 0, 255),
       new StsColor(168, 0, 255, 0, 255),
       new StsColor(254, 0, 0, 255, 255)
 */
    };
    static StsColor[] autoCADSpectrumColors = new StsColor[]
              {
		new StsColor(0,0,0,0, 255),
		new StsColor(1, 255,0,0, 255),
		new StsColor(2, 255,255,0, 255),
		new StsColor(3, 0,255,0, 255),
		new StsColor(4, 0,255,255, 255),
		new StsColor(5, 0,0,255, 255),
		new StsColor(6, 255,0,255, 255),	
		new StsColor(7, 255,255,255, 255),
		new StsColor(8, 65,65,65, 255),
		new StsColor(9, 128,128,128, 255),
		new StsColor(10, 255,0,0, 255),
		new StsColor(11, 255,170,170, 255),
		new StsColor(12, 189,0,0, 255),
		new StsColor(13, 189,126,126, 255),
		new StsColor(14, 129,0,0, 255),
		new StsColor(15, 129,86,86, 255),
		new StsColor(16, 104,0,0, 255),
		new StsColor(17, 104,69,69, 255),
		new StsColor(18, 79,0,0, 255),
		new StsColor(19, 79,53,53, 255),
		new StsColor(20, 255,63,0, 255),
		new StsColor(21, 255,191,170, 255),
		new StsColor(22, 189,46,0, 255),
		new StsColor(23, 189,141,126, 255),
		new StsColor(24, 129,31,0, 255),   	
		new StsColor(25, 129,96,86, 255),
		new StsColor(26, 104,25,0, 255),
		new StsColor(27, 104,78,69, 255),
		new StsColor(28, 79,19,0, 255),
		new StsColor(29, 79, 59, 53, 255),
		new StsColor(30, 255,127,0, 255),
		new StsColor(31, 255,212,170, 255),	
		new StsColor(32, 189,94,0, 255),
		new StsColor(33, 189,157,126, 255),
		new StsColor(34, 129,64,0, 255),
		new StsColor(35, 129,107, 86, 255),
		new StsColor(36, 104, 52, 0, 255),
		new StsColor(37, 104, 86, 69, 255),
		new StsColor(38, 79,39,0, 255),
		new StsColor(39, 79,66,53, 255),
		new StsColor(40, 255,191, 0, 255),
		new StsColor(41, 255,234,170, 255),
		new StsColor(42, 189,141,0, 255),
		new StsColor(43, 189, 173, 126, 255),
		new StsColor(44, 129, 96, 0, 255),
		new StsColor(45, 129, 118, 86, 255),
		new StsColor(46, 104,78, 0, 255),
		new StsColor(47, 104, 95, 69, 255),
		new StsColor(48, 79, 59, 0, 255),
		new StsColor(49, 79, 73, 53, 255),   	
		new StsColor(50, 255, 255, 0, 255),
		new StsColor(51,255,255,170 , 255),
		new StsColor(52, 189,189,0, 255),
		new StsColor(53, 189,189,126, 255),
		new StsColor(54, 129,129,0, 255),
		new StsColor(55, 129,129,86, 255),
		new StsColor(56, 104,104,0, 255),	
		new StsColor(57, 104,104,69, 255),
		new StsColor(58, 79,79,0, 255),
		new StsColor(59, 79,79,53, 255),
		new StsColor(60, 191,255,0, 255),
		new StsColor(61, 234,255,170, 255),
		new StsColor(62, 141,189,0, 255),
		new StsColor(63, 173, 189,0, 255),
		new StsColor(64, 96,129,0, 255),
		new StsColor(65, 118, 129, 86, 255),
		new StsColor(66, 78,104,0, 255),
		new StsColor(67, 95,104,69, 255),
		new StsColor(68, 59,79,0, 255),
		new StsColor(69, 73, 79, 53, 255),
		new StsColor(70, 127, 255, 0, 255),
		new StsColor(71, 212,255,170, 255),
		new StsColor(72, 94, 189, 0, 255),
		new StsColor(73, 157,189,126, 255),
		new StsColor(74, 64,129,0, 255),   	
		new StsColor(75, 107,129,86, 255),
		new StsColor(76, 52,104,0, 255),
		new StsColor(77, 86,104,69, 255),
		new StsColor(78, 39,79,0, 255),
		new StsColor(79, 66, 79,53, 255),
		new StsColor(80, 63, 255,0, 255),
		new StsColor(81, 191, 255, 170, 255),	
		new StsColor(82, 46,189,0, 255),
		new StsColor(83, 141, 189, 126, 255),
		new StsColor(84, 31,129, 0, 255),
		new StsColor(85, 96, 129, 86, 255),
		new StsColor(86, 25, 104, 0, 255),
		new StsColor(87, 78, 104, 69, 255),
		new StsColor(88, 19,79, 0, 255),
		new StsColor(89, 59, 79, 53, 255),
		new StsColor(90, 0, 255, 0, 255),
		new StsColor(91, 170,255,170, 255),
		new StsColor(92, 0, 189,0, 255),
		new StsColor(93, 126, 189, 126, 255),
		new StsColor(94, 0, 129, 0, 255),
		new StsColor(95, 86, 129,86, 255),
		new StsColor(96, 0, 104, 0, 255),
		new StsColor(97, 69,104,69, 255),
		new StsColor(98, 0,79,0, 255),
		new StsColor(99, 53,79,53, 255),   	
		new StsColor(100, 0,255,63, 255),
		new StsColor(101, 170,255,191, 255),
		new StsColor(102, 0,189,46, 255),
		new StsColor(103, 126,189,141, 255),
		new StsColor(104, 0,129,31, 255),
		new StsColor(105, 86,129,96, 255),
		new StsColor(106, 0,104,25, 255),
		new StsColor(107, 69,104,78, 255),
		new StsColor(108, 0,79,19, 255),
		new StsColor(109, 53,79,59, 255),
		new StsColor(110, 0,255,127, 255),
		new StsColor(111, 170,255,212, 255),
		new StsColor(112, 0,189,94, 255),
		new StsColor(113, 126,189,157, 255),
		new StsColor(114, 0,129,64, 255),
		new StsColor(115, 86,129,107, 255),
		new StsColor(116, 0,104,52, 255),
		new StsColor(117, 69,104,86, 255),
		new StsColor(118, 0,79,39, 255),
		new StsColor(119, 53,79,66, 255),
		new StsColor(120, 0,255,191, 255),
		new StsColor(121, 170,255,191, 255),
		new StsColor(122, 0,189,141, 255),
		new StsColor(123, 126,189,173, 255),
		new StsColor(124, 0,129,96, 255),
		new StsColor(125, 86,129,118, 255),
    			new StsColor(126, 0,104,78, 255),
    			new StsColor(127, 69,104,95, 255),
    			new StsColor(128, 0,79,59, 255),
    			new StsColor(129, 53,79,73, 255),
    			new StsColor(130, 0,255,255, 255),
    			new StsColor(131, 170,255,255, 255),	
    			new StsColor(132, 0,189,189, 255),
    			new StsColor(133, 126,189,189, 255),
    			new StsColor(134, 0,129,129, 255),
    			new StsColor(135, 86,129,129, 255),
    			new StsColor(136, 0,104,104, 255),
    			new StsColor(137, 69,104,104, 255),
    			new StsColor(138, 0,79,79, 255),
    			new StsColor(139, 53,79,79, 255),
    			new StsColor(140, 0,191,255, 255),
    			new StsColor(141, 170,234,255, 255),
    			new StsColor(142, 0,141,189, 255),
    			new StsColor(143, 126,173,189, 255),
    			new StsColor(144, 0,96,129, 255),
    			new StsColor(145, 86,118,129, 255),
    			new StsColor(146, 0,78,104, 255),
    			new StsColor(147, 69,95,104, 255),
				new StsColor(148, 0,59,79, 255),
				new StsColor(149, 53,73, 79, 255),   	
    			new StsColor(150, 0,127,255, 255),
    			new StsColor(151, 170,212,255, 255),
    			new StsColor(152, 0,94,189, 255),
    			new StsColor(153, 126,157,189, 255),
    			new StsColor(154, 0,64,129, 255),
    			new StsColor(155, 86,107,129, 255),
    			new StsColor(156, 0,52,104, 255),
    			new StsColor(157, 69,86,104, 255),
    			new StsColor(158, 0,39,79, 255),
    			new StsColor(159,53,66,79 , 255),
    			new StsColor(160, 0,63,255, 255),
    			new StsColor(161, 170,191,255, 255),
    			new StsColor(162, 0,46,189, 255),
    			new StsColor(163, 126,141,189, 255),
    			new StsColor(164, 0,31,129, 255),
    			new StsColor(165, 86,96,129, 255),
    			new StsColor(166, 0,25,104, 255),
    			new StsColor(167, 69,78,104, 255),
    			new StsColor(168, 0,19,79, 255),
    			new StsColor(169, 53,59,79, 255),
    			new StsColor(170, 0,0,255, 255),
    			new StsColor(171, 170,170,255, 255),
    			new StsColor(172, 0,0,189, 255),
    			new StsColor(173, 126,126,189, 255),
    			new StsColor(174, 0,0,129, 255),
    			new StsColor(175, 86, 86, 129, 255),
    			new StsColor(176, 0, 0, 104, 255),
    			new StsColor(177, 69, 69, 104, 255),
    			new StsColor(178, 0, 0, 79, 255),
    			new StsColor(179, 53, 53, 79, 255),
    			new StsColor(181, 191, 170, 255, 255),
    			new StsColor(182, 46, 0, 189, 255),
    			new StsColor(183, 141, 126, 189, 255),
    			new StsColor(184, 31, 0, 129, 255),
    			new StsColor(185, 96, 86, 129, 255),
    			new StsColor(186, 25, 0, 104, 255),
    			new StsColor(187, 78, 69, 104, 255),
    			new StsColor(188, 19, 0, 79, 255),
    			new StsColor(189, 59, 53, 79, 255),
    			new StsColor(190, 127, 0, 255, 255),
    			new StsColor(191, 212, 170, 255, 255),
    			new StsColor(192, 94, 0, 189, 255),
    			new StsColor(193, 157, 126, 189, 255),
    			new StsColor(194, 64, 0, 129, 255),
    			new StsColor(195, 107, 86, 129, 255),
    			new StsColor(196, 52, 0, 104, 255),
    			new StsColor(197, 86, 69, 104, 255),
    			new StsColor(198, 39, 0, 79, 255),
    			new StsColor(199, 66, 53, 79, 255),
    			new StsColor(200, 191, 0, 255, 255),
    			new StsColor(201, 234, 170, 255, 255),
    			new StsColor(202, 141, 0, 189, 255),
    			new StsColor(203, 173, 126, 189, 255),
    			new StsColor(204, 96, 0, 129, 255),
    			new StsColor(205, 118, 86, 129, 255),
    			new StsColor(206, 78, 0, 104, 255),
    			new StsColor(207, 95, 69, 104, 255),
    			new StsColor(208, 59,  0, 79, 255),
    			new StsColor(209, 73, 53, 79, 255),
    			new StsColor(210, 255, 0, 255, 255),
    			new StsColor(211, 255, 170, 255, 255),
    			new StsColor(212, 189, 0, 189, 255),
    			new StsColor(213, 189, 126, 189, 255),
    			new StsColor(214, 129, 0, 129, 255),
    			new StsColor(215, 129, 86, 129, 255),
    			new StsColor(216, 104, 0, 104, 255),
    			new StsColor(217, 104, 69, 104, 255),
    			new StsColor(218, 79, 0, 79, 255),
    			new StsColor(219, 79, 53, 79, 255),
    			new StsColor(220, 255, 0, 191, 255),
    			new StsColor(221, 255, 170, 234, 255),
    			new StsColor(222, 189, 0, 141, 255),
    			new StsColor(223, 189, 126, 173, 255),
    			new StsColor(224, 129, 0, 96, 255),
    			new StsColor(225, 129, 86, 118, 255),
    			new StsColor(226, 104, 0, 78, 255),
    			new StsColor(227, 104, 69, 95, 255),
    			new StsColor(228, 79, 0, 59, 255),
    			new StsColor(229, 79, 53, 73, 255),
    			new StsColor(230, 255,0,127, 255),
    			new StsColor(231, 255, 0, 127, 255),
    			new StsColor(232, 189, 0, 94, 255),
    			new StsColor(233, 189, 126, 157, 255),
    			new StsColor(234, 129, 0, 64, 255),
    			new StsColor(235, 129, 86, 107, 255),
    			new StsColor(236, 104, 0, 52, 255),
    			new StsColor(237, 104, 69, 86, 255),
    			new StsColor(238, 79, 0, 39, 255),
    			new StsColor(239, 79, 53, 66, 255),
    			new StsColor(240, 255, 0, 63, 255),
    			new StsColor(241, 255, 170, 191, 255),
    			new StsColor(242, 189, 0, 46, 255),
    			new StsColor(243, 189, 126, 141, 255),
    			new StsColor(244, 129, 0, 31, 255),
    			new StsColor(245, 129, 86, 96, 255),
    			new StsColor(246, 104, 0, 25, 255),
    			new StsColor(247, 104, 69, 78, 255),
    			new StsColor(248, 79, 0, 19, 255),
    			new StsColor(249, 79, 53,59, 255),
    			new StsColor(250, 51,51,51, 255),
    			new StsColor(251, 80, 80, 80, 255),
    			new StsColor(252, 105,105,105, 255),
    			new StsColor(253, 130,130,130, 255),
    			new StsColor(254, 190, 190, 190, 255),
    			new StsColor(255, 255, 255, 255, 255),    			
              };

    static public final String SPECTRUM_RWB = "RedWhiteBlue";
    static public final String SPECTRUM_BWYR = "BlackWhiteYellowRed";
    static public final String SPECTRUM_SIMULATION = "Simulation Spectrum";
    static public final String SPECTRUM_PROPERTIES = "Properties Spectrum";
    static public final String SPECTRUM_RAINBOW = "Rainbow Spectrum";
    static public final String SPECTRUM_SEMBLANCE = "Semblance Spectrum";
    static public final String SPECTRUM_GRAYSCALE = "Grayscale Spectrum";
    static public final String SPECTRUM_AUTOCAD = "AutoCAD";
    static public final String SPECTRUM_STACKS = "Default Stacks Spectrum";
    static public final String SPECTRUM_BWY = "BlueWhiteYellow";
    static public final String SPECTRUM_TEMPERATURES = "Default Temperature Colors";
    static public final String SPECTRUM_INTERVAL_VEL = "Default Interval Velocity Spectrum";
    static public final String SPECTRUM_RMS_VEL = "Default RMS Velocity Spectrum";
    static public final String SPECTRUM_WATER_SAT = "Default Water Sat Spectrum";

    static public String[] cannedSpectrums =
    {
       SPECTRUM_RWB, SPECTRUM_RAINBOW, SPECTRUM_BWYR, SPECTRUM_SIMULATION,
       SPECTRUM_PROPERTIES, SPECTRUM_SEMBLANCE, SPECTRUM_GRAYSCALE, SPECTRUM_AUTOCAD
    };

    static public String[] velocitySpectrums = { SPECTRUM_INTERVAL_VEL, SPECTRUM_RMS_VEL,  SPECTRUM_INTERVAL_VEL, SPECTRUM_INTERVAL_VEL};

    //TODO create spectrumDefinitions: class StsSpectrumDefinition(name, colors).
    //TODO For example, StsSpectrumDefinition temperature = new StsSpectrumDefinition(SPECTRUM_TEMPERATURES, defaultTempColors)
    //TODO create temperature spectrums in this class, they are made persistent when "adopted" by some parent object, who might subsequently edit them,
    //TODO so they need to be copies to a persisted instance.
    public StsSpectrumClass()
    {
        userName = "Color Spectrum";
    }

    /** Called to create a default set of spectrums */
    public boolean createSpectrums()
    {
        try
        {
            StsSpectrum spectrum = new StsSpectrum("Basic");
            spectrum.setBasic32Colors();
            spectrum.addToModel();

            spectrum = new StsSpectrum(SPECTRUM_RWB);
            StsColor[] colors = {StsColor.BLUE, StsColor.WHITE, StsColor.RED};
            spectrum.setInterpColors(colors, 255);
            spectrum.addToModel();

            spectrum = StsSpectrum.createRainbowSpectrum(SPECTRUM_RAINBOW, 255);
            spectrum.addToModel();

            spectrum = new StsSpectrum("RGB10");
            spectrum.setBasic10Colors();
            spectrum.addToModel();

            spectrum = new StsSpectrum(StsSpectrumClass.SPECTRUM_BWYR);
            spectrum.setInterpColorsWithIdx(defaultSeismicSpectrumColors, 255);
            spectrum.addToModel();

            spectrum = new StsSpectrum(SPECTRUM_SIMULATION);
            spectrum.setInterpColorsWithIdx(defaultSimulationSpectrumColors, 255);
            spectrum.addToModel();

            spectrum = new StsSpectrum(SPECTRUM_PROPERTIES);
            spectrum.setInterpColorsWithIdx(defaultPropertySpectrumColors, 255);
            spectrum.addToModel();

            spectrum = new StsSpectrum(SPECTRUM_SEMBLANCE);
            spectrum.setInterpColorsWithIdx(defaultSemblanceColors, 255);
            spectrum.addToModel();
            
            spectrum = new StsSpectrum(SPECTRUM_INTERVAL_VEL);
            spectrum.setInterpColorsWithIdx(defaultSemblanceColors, 255);
            spectrum.addToModel();            

            spectrum = new StsSpectrum(SPECTRUM_RMS_VEL);
            spectrum.setInterpColorsWithIdx(defaultSemblanceColors, 255);
            spectrum.addBlackContours(25);
            spectrum.addToModel();

            spectrum = new StsSpectrum(SPECTRUM_WATER_SAT);
            spectrum.setInterpColorsWithIdx(defaultWSatColors, 255);
            spectrum.addToModel();
            
            spectrum = new StsSpectrum(SPECTRUM_STACKS);
            spectrum.setInterpColorsWithIdx(defaultSeismicSpectrumColors, 255);
            spectrum.addToModel();

            spectrum = new StsSpectrum(SPECTRUM_GRAYSCALE);
            spectrum.setGrayScale(255);
            spectrum.addToModel();

            spectrum = new StsSpectrum(SPECTRUM_AUTOCAD);
            spectrum.setStsColors(autoCADSpectrumColors);
            spectrum.addToModel();

            spectrum = new StsSpectrum(SPECTRUM_TEMPERATURES);
            spectrum.setInterpColorsWithIdx(defaultTempColors, 255);
            spectrum.addToModel();
            
            spectrum = new StsSpectrum(SPECTRUM_BWY);  //colorscale based of Focus blue-white-yellow
            StsColor[] bwy = {new StsColor(255, 229, 5), StsColor.WHITE, new StsColor(2, 8, 78)};
            spectrum.setInterpColors(bwy, 255);
            spectrum.addToModel();
            
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSpectrum.createSpectrums() failed.",
                    e, StsException.WARNING);
            return false;
        }
    }

    public StsSpectrum getSpectrum(String name) { return (StsSpectrum)getObjectWithName(name); }

    public void incrementSpectrumColor(String name)
    {
        StsSpectrum spectrum = (StsSpectrum)getObjectWithName(name);
        spectrum.incrementCurrentColor();
    }

    public StsColor getCurrentSpectrumColor(String name)
    {
        StsSpectrum spectrum = (StsSpectrum)getObjectWithName(name);
        if(spectrum == null) return StsColor.RED;
        return spectrum.getCurrentColor();
    }
}
