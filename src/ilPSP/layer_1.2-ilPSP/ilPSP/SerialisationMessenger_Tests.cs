using System;
using MessagePack;
using NUnit.Framework;

namespace ilPSP.Utils {

/*
    /// <summary>Compatibility checks for the MPI object serializer.</summary>
    [TestFixture]
    public static class SerialisationMessengerTests {
        /// <summary>
        /// Agglomeration messages contain vectors inside arrays and value tuples;
        /// default vectors, dimensions and floating-point bits must be retained.
        /// </summary>
        [Test]
        public static void VectorPayloadsPreserveState() {
            var vectors = new[] {
                default(Vector), new Vector(-0.0), new Vector(1.25, double.PositiveInfinity),
                new Vector(double.NaN, double.NegativeInfinity, Math.PI)
            };
            var input = new VectorMessage {
                Vectors = vectors,
                Pair = Tuple.Create(1234567890123L, vectors[2]),
                Boxed = vectors[3],
                Nullable = vectors[1]
            };
            var bytes = MessagePackSerializer.Serialize<object>(input, SerialisationMessenger.SerializationOptions);
            var output = MessagePackSerializer.Deserialize<object>(bytes, SerialisationMessenger.SerializationOptions) as VectorMessage;
            Assert.IsNotNull(output);
            Assert.AreEqual(input.Vectors.Length, output.Vectors.Length);
            for (int i = 0; i < vectors.Length; i++) Check(vectors[i], output.Vectors[i]);
            Assert.AreEqual(input.Pair.Item1, output.Pair.Item1);
            Check(input.Pair.Item2, output.Pair.Item2);
            Assert.IsInstanceOf<Vector>(output.Boxed);
            Check((Vector)input.Boxed, (Vector)output.Boxed);
            Assert.IsTrue(output.Nullable.HasValue);
            Check(input.Nullable.Value, output.Nullable.Value);
        }

        static void Check(Vector expected, Vector actual) {
            Assert.AreEqual(expected.Dim, actual.Dim);
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected.x), BitConverter.DoubleToInt64Bits(actual.x));
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected.y), BitConverter.DoubleToInt64Bits(actual.y));
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected.z), BitConverter.DoubleToInt64Bits(actual.z));
            Assert.AreEqual(0, actual.Dummy_256bitAlign);
        }

        sealed class VectorMessage {
            public Vector[] Vectors;
            public Tuple<long, Vector> Pair;
            public object Boxed;
            public Vector? Nullable;
        }
    }
*/
}
