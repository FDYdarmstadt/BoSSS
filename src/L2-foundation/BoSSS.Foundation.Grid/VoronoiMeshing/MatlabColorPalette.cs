using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    struct RGBColor
    {
        public byte Red;
        public byte Green;
        public byte Blue;
    }

    static class ColorPalette
    {
        public static RGBColor Silver()
        {
            return new RGBColor
            {
                Red = 192,
                Green = 192,
                Blue = 192
            };
        }

        public static RGBColor Slategrey()
        {
            return new RGBColor
            {
                Red = 112,
                Green = 128,
                Blue = 144
            };
        }

        public static RGBColor Lightblue()
        {
            return new RGBColor
            {
                Red = 102,
                Green = 153,
                Blue = 255
            };
        }
    }

    static class MatlabColorPalette
    {
        public static MatlabColor Silver()
        {
            return new MatlabColor(ColorPalette.Silver());
        }

        public static MatlabColor Slategrey()
        {
            return new MatlabColor(ColorPalette.Slategrey());
        }

        public static MatlabColor Lightblue()
        {
            return new MatlabColor(ColorPalette.Lightblue());
        }
    }

    class MatlabColor
    {
        public string color;

        public CultureInfo cultureInfo = CultureInfo.InvariantCulture;

        public MatlabColor() { }

        public MatlabColor(RGBColor color)
        {
            string red = Normalized(color.Red).ToString(cultureInfo);
            string green = Normalized(color.Green).ToString(cultureInfo); ;
            string blue = Normalized(color.Blue).ToString(cultureInfo); ;
            this.color = $"[{red} {green} {blue}]";
        }

        double Normalized(byte color)
        {
            return color / 255.0;
        }

        public static implicit operator string(MatlabColor color)
        {
            return color.color;
        }
    }
}
