using System;
using System.Globalization;
using Newtonsoft.Json;

namespace ilPSP.Utils {
    /// <summary>
    /// Retains date types in object-valued MPI payloads while leaving date-like
    /// strings unchanged. JSON has no distinct date token, unlike BSON.
    /// </summary>
    internal sealed class TypePreservingDateConverter : JsonConverter {
        public override bool CanRead => false;

        public override bool CanConvert(Type objectType) {
            return objectType == typeof(DateTime) || objectType == typeof(DateTimeOffset);
        }

        public override void WriteJson(JsonWriter writer, object value, JsonSerializer serializer) {
            writer.WriteStartObject();
            writer.WritePropertyName("$type");
            writer.WriteValue(value.GetType().AssemblyQualifiedName);
            writer.WritePropertyName("$value");
            writer.WriteValue(((IFormattable)value).ToString("O", CultureInfo.InvariantCulture));
            writer.WriteEndObject();
        }

        public override object ReadJson(JsonReader reader, Type objectType, object existingValue, JsonSerializer serializer) {
            // Json.NET reads the $type/$value envelope with its built-in primitive reader.
            throw new NotSupportedException();
        }
    }
}
