using BlazorApp1.Shared.Enums;

namespace BlazorApp1.Shared.Models;

public record TaskModel(Guid Id, string Title, string Description, Status Status);
