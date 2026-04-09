package Bio::Graphics;

use strict;
use warnings;

package Bio::Graphics::Panel;

sub new {
    my ($class, %args) = @_;
    return bless {
        args => \%args,
        tracks => [],
    }, $class;
}

sub add_track {
    my ($self, @args) = @_;
    my $track = bless {
        args => \@args,
        features => [],
    }, 'Bio::Graphics::Track';
    push @{$self->{tracks}}, $track;
    return $track;
}

sub image_and_map {
    my ($self, %args) = @_;
    my $mapname = $args{-url} || 'virusdetect-map';
    return ('about:blank', '', $mapname);
}

package Bio::Graphics::Track;

sub add_feature {
    my ($self, $feature) = @_;
    push @{$self->{features}}, $feature;
    return 1;
}

package Bio::Graphics::Feature;

sub new {
    my ($class, %args) = @_;
    return bless {
        args => \%args,
        tags => {},
    }, $class;
}

sub has_tag {
    my ($self, $tag) = @_;
    return exists $self->{tags}{$tag};
}

sub each_tag_value {
    my ($self, $tag) = @_;
    return () unless exists $self->{tags}{$tag};
    my $value = $self->{tags}{$tag};
    return ref($value) eq 'ARRAY' ? @$value : ($value);
}

1;
